#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include <chrono>
#include <pthread.h>
#include <omp.h>
#include <fstream>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

double binomialCoefficients(int n, int k) {
   double C[k+1];
   memset(C, 0, sizeof(C));
   C[0] = 1;
   for (int i = 1; i <= n; i++) {
      for (int j = std::min(i, k); j > 0; j--)
         C[j] = C[j] + C[j-1];
   }
   return C[k];
}

class Quasi1d {       // The class
  public:             // Access specifier
    int nsize;
    double eps;
    double bp;
    double amin;

    double Lmean;
    Eigen::Matrix<double, Eigen::Dynamic, 1> xvalues;




    Quasi1d(int nsize0, double eps0){
        nsize = nsize0;
        eps = eps0;
        amin = sqrt(1-eps*eps);

        //Internal stuff
        am.resize(nsize, nsize);
        M.resize(nsize, nsize);
        Qt.resize(nsize,nsize);
        Qtm.resize(nsize,nsize);
        Qt2.resize(nsize,nsize);
        Qtm2.resize(nsize,nsize);
        UpdateMatrixA();
        
    }

    int SetPressure(double bp0){
        //Create matrix M
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                M(i,j) = omega(bp0,i,j);
            }
        }

        double eigval;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es;
        es.compute(M);
        xvalues = es.eigenvectors().col(nsize-1);
        Lmean = 1.0/es.eigenvalues()(nsize-1);

        //Prepare values
        double sum = 0;
        for(int i = 0; i<nsize; i++){
            xvalues(i) = xvalues(i)*xvalues(i);
            sum += xvalues(i);
        }
        xvalues /= sum;

        //Update internal bp
        std::cout << "Updating internal pressure value." << std::endl;
        bp = bp0;
        return 0;

    }

    // Compute density of the system at a certain bp
    double Density(double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double nsum = 0;
        //Igual puedo optimizar este bucle y hacer la mitad del bucle
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                nsum += xvalues(i)*xvalues(j)*Lmean*omegap(bp,i,j)/sqrt(xvalues(i)*xvalues(j));
            }
        }
        return -1.0/nsum;
    }

    // Compute Gir(r) analiticamente para un nmax
    double Gijr(double x, double bp0, int nmax, int i, int j){
        if (bp0 != bp)  SetPressure(bp0);
        //Optimizacion: vemos en que capa de proximos vecinos cae x y ajustamos el valor de nmax si es necesario
        int nshell = int(x/am(i,j));
        if (nshell == 0) return 0;

        if (nshell < nmax){
            nmax = nshell; //keep lowest n value because they are going to yield the same result
            //std::cout << "Using nshell = " << nmax << std::endl; 
        }
        double gij = 0.0;
        for (int ni=0; ni<nmax; ni++){
            gij += pow(Lmean, ni+1) * qij(x,bp,ni+1,i,j);
        }
        gij /=(Density(bp0)*sqrt(xvalues(i)*xvalues(j)));
        return gij;
    }

    // Compute global G(r)
    double Gr(double x, double bp0, int nmax){
        if (bp0 != bp)  SetPressure(bp0);
        //Optimizacion: vemos en que capa de proximos vecinos cae x y ajustamos el valor de nmax si es necesario
        int nshell = int(x/amin);
        if (nshell == 0) return 0;

        if (nshell < nmax){
            nmax = nshell; //keep lowest n value because they are going to yield the same result
            //std::cout << "Using nshell = " << nmax << std::endl; 
        }
        double gtotal = 0.0;
        #pragma omp parallel for reduction(+:gtotal) num_threads(16) schedule(dynamic)
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                if (ix!=jx){
                    gtotal += 2.0*xvalues(ix)*xvalues(jx)*Gijr(x,bp0,nmax,ix,jx);
                }
                else{
                    gtotal = gtotal + xvalues(ix)*xvalues(jx)*Gijr(x,bp0,nmax,ix,jx);
                } 
            }

        }
        return gtotal;
    }

    // Compute Gs real numbers
    double Gs(double s, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double gtotal = 0;
        //Create matrix Qt
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                Qt(i,j) = Lmean*omega(s+bp0,i,j);
                Qtm(i,j) = -Qt(i,j);
                if(i==j) Qtm(i,j) = 1.0+Qtm(i,j);

            }
        }
        //Invert matrix and all that
        Qt = Qt*(Qtm.inverse());

        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                if (i!=j){
                    gtotal += 2.0*sqrt(xvalues(i)*xvalues(j))*Qt(i,j);
                }
                else{
                    gtotal += sqrt(xvalues(i)*xvalues(j))*Qt(i,j);
                } 
            }
        }
        return gtotal/Density(bp);
    }

    // Compute Gs with complex numbers, I need to change that, of course
    std::complex<double> Gs2(std::complex<double> s,  double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        std::complex<double> gtotal = 0;

        //Create matrix Qt
        #pragma omp parallel for num_threads(16)
        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                Qt2(i,j) = Lmean*omega2(s+bp0,i,j);
                Qtm2(i,j) = -Qt2(i,j);
                if(i==j) Qtm2(i,j) = std::complex<double>(1.0,0.0)+Qtm2(i,j);

                // For symmetric indices too
                Qt2(j,i) = Qt2(i,j);
                Qtm2(j,i) = Qtm2(i,j);

            }
        }
        //Invert matrix and all that
        Qt2 = Qt2*(Qtm2.inverse());

        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                if (i!=j){
                    gtotal += 2.0*sqrt(xvalues(i)*xvalues(j))*Qt2(i,j);
                }
                else{
                    gtotal += sqrt(xvalues(i)*xvalues(j))*Qt2(i,j);
                } 
            }
        }
        std::complex<double> density(Density(bp),0);
        return gtotal/density;
    }

    // Compute Gs with complex numbers, I need to change that, of course
    std::complex<double> Gijs2(std::complex<double> s,  double bp0, int ii, int jj){
        if (bp0 != bp)  SetPressure(bp0);
        std::complex<double> gtotal = 0;

        //Create matrix Qt
        #pragma omp parallel for num_threads(16)
        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                Qt2(i,j) = Lmean*omega2(s+bp0,i,j);
                Qtm2(i,j) = -Qt2(i,j);
                if(i==j) Qtm2(i,j) = std::complex<double>(1.0,0.0)+Qtm2(i,j);

                // For symmetric indices too
                Qt2(j,i) = Qt2(i,j);
                Qtm2(j,i) = Qtm2(i,j);

            }
        }
        //Invert matrix and all that
        Qt2 = Qt2*(Qtm2.inverse());

        //for (int i=0; i<nsize; i++){
        //    for(int j=i; j<nsize; j++){
                //if (i!=j){
                //    gtotal += 2.0*sqrt(xvalues(i)*xvalues(j))*Qt2(i,j);
                //}
                //else{
        gtotal += Qt2(ii,jj)/sqrt(xvalues(ii)*xvalues(jj));
                //} 
            //}
        //}
        std::complex<double> density(Density(bp),0);
        return gtotal/density;
    }

    //Compute Sk
    double Sk(double k, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double stotal = 0;

        return 1.0 + Density(bp0)*std::real(Gs2(std::complex<double>(0.0,k),bp0)+Gs2(std::complex<double>(0.0,-k),bp0));

    }

    //Gr calculated from numerical Laplace inverse of G(s)
    double GrInv(double t, double bp0, int ntr, int meuler){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::complex<double> suma = Gs2(xx,bp0)/2.0;

        std::complex<double> argument(0,0);
        for (int i=1;i<ntr+1;i++){
            suma = suma + std::complex<double>(pow(-1,i),0) * Gs2(std::complex<double>(xx, i*hh),bp0);
        }

        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        su[0] = suma;

        //#avgsu
        std::complex<double> avgsu(0,0);
        std::complex<double> avgsu1(0,0);
        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            su[i+1] = su[i] + pow(-1,nn) * Gs2(std::complex<double>(xx,nn*hh),bp0);
            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            avgsu += bn*su[i];
            avgsu1 += bn*su[i+1];

        }
        //for(int i=0;i<meuler+1;i++){
        //    avgsu += binomialCoefficients(meuler,(i+1)-1)*su[i];
        //    avgsu1 += binomialCoefficients(meuler,(i+1)-1)*su[i+1];
        //    //#print(i, avgsu)
        //}
        double pp = pow(2.0,meuler);
        std::complex<double> fun = uu*avgsu/(pp);
        std::complex<double> fun1 = uu*avgsu1/(pp);
        return real(fun1);
        //return 0.0;

    }

        //Gr calculated from numerical Laplace inverse of G(s)
    double GijrInv(double t, double bp0, int ntr, int meuler, int ii, int jj){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::complex<double> suma = Gijs2(xx,bp0, ii, jj)/2.0;

        std::complex<double> argument(0,0);
        for (int i=1;i<ntr+1;i++){
            suma = suma + std::complex<double>(pow(-1,i),0) * Gijs2(std::complex<double>(xx, i*hh),bp0, ii, jj);
        }

        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        su[0] = suma;

        //#avgsu
        std::complex<double> avgsu(0,0);
        std::complex<double> avgsu1(0,0);
        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            su[i+1] = su[i] + pow(-1,nn) * Gijs2(std::complex<double>(xx,nn*hh),bp0, ii, jj);
            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            avgsu += bn*su[i];
            avgsu1 += bn*su[i+1];

        }
        double pp = pow(2.0,meuler);
        std::complex<double> fun = uu*avgsu/(pp);
        std::complex<double> fun1 = uu*avgsu1/(pp);
        return real(fun1);
        //return 0.0;

    }

    double ProbDensity(double x, int n, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        
        double gtotal = 0.0;
        #pragma omp parallel for reduction(+:gtotal) num_threads(16) schedule(dynamic)
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                if (ix!=jx){
                    gtotal += 2.0*sqrt(xvalues(jx)*xvalues(ix))*pow(Lmean,n)*qij(x, bp0, n, ix, jx);
                }
                else{
                    gtotal += sqrt(xvalues(jx)*xvalues(ix))*pow(Lmean,n)*qij(x, bp0, n, ix, jx);
                } 
                //std::cout << ix << ", " << jx << ", " << sqrt(xvalues(jx)/xvalues(ix))<< ", " <<Lmean<< ", " <<qij(x, bp0, n, ix, jx)<< ", " << gtotal << std::endl;
            }

        }
        return gtotal;
    }
    
    
   private:

    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> am;
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> M;
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> Qt;
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> Qtm;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qt2;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtm2;

    double omega(double s, int i, int j){
        return exp(-am(i,j)*s)/s;
    }

    std::complex<double> omega2(std::complex<double> s, int i, int j){
        return exp(-am(i,j)*s)/s;
    }

    double omegap(double s, int i, int j){
        return -exp(-am(i,j)*s)*(1.0+am(i,j)*s)/(s*s);
    }

    double y(int i){
        return -eps/2.0 + ((double)i-1)*eps/((double)nsize-1);
    }

    int UpdateMatrixA(){
        double temp;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                temp = (y(i)-y(j))*(y(i)-y(j));
                am(i,j) = sqrt(1.0-temp);
            }
        }
        return 0;
    }  

    double qij(double x, double bp0, int n, int ii, int jj){
      if (bp0 != bp)  SetPressure(bp0);
        if (n>1){

            std::vector<int> kvec(n+1,0); 
            kvec[0]=ii;
            kvec[n]=jj;

            int p = 1; //Used to increment all of the indicies correctly, at the end of each loop.
            double stotal = 0.0;
            while (kvec[n]==jj) {

                stotal += qikj(x,bp0,kvec,n);
                //Update indices correctly
                kvec[1]++;
                while(kvec[p]==nsize && p<n) {
                    kvec[p]=0;
                    kvec[++p]++; //increase p by 1, and increase the next (p+1)th index
                    if(kvec.at(p)!=nsize)
                        p=1;
                }
            }
         return stotal;
      }
      else {
        std::vector<int> kvec{ii,jj};
        return qikj(x,bp0,kvec,1);
      }

    }

    double qikj(double x, double bp0, std::vector<int> kvec, int n){
        if (bp0 != bp)  SetPressure(bp0);
        if (kvec.size() != n+1){
            std::cout << "Wrong size of kvec! Something went wrong with the input parameters \n";
            return 0.0;
        }
        double sa = x;
        for (int h=0; h<kvec.size()-1; h++){
            sa -= am(kvec.at(h),kvec.at(h+1));
        }
        if (sa > 0){
            return exp(-bp0*x)*pow(sa, n-1.0)/tgamma(n); //tgamma(n+1)=factorial(n)
        }
        else{
            return 0.0;
        }
    }
};



int main(int argc, char* argv[]) {


    double eps = atof(argv[1]);
    int nsize = atoi(argv[3]);
    int nmax = 4;
    Quasi1d f(nsize, eps);  // Create an object of MyClass

    double bp = atof(argv[2]);
    f.SetPressure(bp);

    //Compute Sk
    if(0){
        double rmin = 4;
        double rmax = 8;
        int rcount = 300;
        double rstep = (rmax-rmin)/rcount;

        double ri;
        std::vector<double> rvalues(rcount,0.0);
        std::vector<double> gvalues(rcount,0.0);

  
        auto begin = std::chrono::high_resolution_clock::now();
        int bwidth = 70;
        double progress;
        double skmax = 0;
        double pkmax = 0;
        for(int k=0;k<=rcount;k++){
            ri = rmin + k*rstep;
            rvalues[k] = ri;
            gvalues[k] = f.Sk(ri, bp);
            //std::cout << gvalues[k] << ", " << skmax << std::endl;
            if (gvalues[k] > skmax){
                //std::cout << "Updating SkMax" << std::endl;
                skmax = gvalues[k];
                pkmax = ri;
            }

            // //Progress bar
            // progress = double(k)/double(rcount);
            // int pos = bwidth*progress;
            // std::cout << "[";
            // for (int i = 0; i < bwidth; ++i) {
            //     if (i < pos) std::cout << "=";
            //     else if (i == pos) std::cout << ">";
            //     else std::cout << " ";
            // }
            // std::cout << "] " << int(progress * 100.0) << " %\r";
            // std::cout.flush();

        }
        std::cout << std::endl;

        std::cout << std::fixed << std::setprecision(5) << "Position of peak: " << f.Density(bp) << " " << bp << " " << pkmax << std::endl;
    
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

        std::cout << "Elapsed time: " << elapsed.count()*1e-9 << std::endl;

        //Write results to file
        std::ofstream ofile;
        ofile.open("Sk_output.txt");
        for(int k=0;k<=rcount;k++){
            ofile << rvalues[k] << "\t" << gvalues[k] << "\n";
    }
    ofile.close();
    }
    else{    //Compute Gr
        double rmin = f.amin;
        double rmax = 5;
        int rcount = 300;
        double rstep = (rmax-rmin)/rcount;

        double ri;
        std::vector<double> rvalues(rcount,0.0);
        std::vector<double> gvalues(rcount,0.0);
        std::vector<double> gvalues11(rcount,0.0);
        std::vector<double> gvalues12(rcount,0.0);
        std::vector<double> gvalues22(rcount,0.0);
        std::vector<double> gvalues13(rcount,0.0);
    
        auto begin = std::chrono::high_resolution_clock::now();
        int bwidth = 70;
        double progress;
        for(int k=0;k<=rcount;k++){
            ri = rmin + k*rstep;
            rvalues[k] = ri;
            //g(r)
            if (ri < (nmax+1)*f.amin){
                //gvalues[k] = f.Gr(ri, bp, nmax);
                gvalues11[k] = f.Gijr(ri, bp, nmax,0,0);
                gvalues13[k] = f.Gijr(ri, bp, nmax,0,nsize-1);
                gvalues12[k] = f.Gijr(ri, bp, nmax,0,int((nsize+1)/2));
                gvalues22[k] = f.Gijr(ri, bp, nmax,int((nsize+1)/2),int((nsize+1)/2));
            }
            else{
                //gvalues[k] = f.GrInv(ri, bp, 100,30);
                gvalues11[k] = f.GijrInv(ri, bp, 100,30,0,0);
                gvalues13[k] = f.GijrInv(ri, bp, 100,30,0,nsize-1);
                gvalues12[k] = f.GijrInv(ri, bp, 100,30,0,int((nsize+1)/2));
                gvalues22[k] = f.GijrInv(ri, bp, 100,30,int((nsize+1)/2),int((nsize+1)/2));
            }

            //Probability Density
            //gvalues[k] = f.ProbDensity(ri,3,bp);

            //Progress bar
            progress = double(k)/double(rcount);
            int pos = bwidth*progress;
            std::cout << "[";
            for (int i = 0; i < bwidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            //std::cout << k << std::endl;
            std::cout.flush();

        }
        std::cout << std::endl;
    
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

        std::cout << "Elapsed time: " << elapsed.count()*1e-9 << std::endl;

        //Write results to file
        std::ofstream ofile;
        ofile.open("output.txt");
        for(int k=0;k<=(int)rcount;k++){
            ofile << rvalues[k] << "\t" << gvalues[k] << " " << gvalues11[k] << " " << gvalues12[k]<< " " << gvalues22[k]<< " " << gvalues13[k] << "\n";
        }
        ofile.close();
    }

  return 0;
}