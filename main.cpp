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
        UpdateMatrixA();

        //Useful matrices for later on
        Qt.resize(nsize,nsize);
        Qtm.resize(nsize,nsize);
        
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

        if (int(x/am(i,j)) == 0)  return 0;

        double gij = 0.0;
        for (int ni=0; ni<nmax; ni++){
            gij += pow(Lmean, ni+1) * qij(x,bp,ni+1,i,j);
        }
        gij /=(Density(bp0)*sqrt(xvalues(i)*xvalues(j)));
        return gij;
    }

    // Compute global G(r)
    std::vector<double> Gr(double x, double bp0, int nmax){
        if (bp0 != bp)  SetPressure(bp0);
        std::vector<double> gtotal = std::vector<double>(5,0.0);

        //Optimizacion: vemos en que capa de proximos vecinos cae x y ajustamos el valor de nmax si es necesario
        int nshell = int(x/amin);
        if (nshell == 0) return gtotal;
        else if (nshell < nmax)  nmax = nshell; //keep lowest n value because they are going to yield the same result

        double gtemp = 0.0;
        double gmean = 0.0;
        
        #pragma omp parallel for reduction(+:gmean) num_threads(16) schedule(dynamic)
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                //gtemp = Gijr(x,bp0,nmax,ix,jx);
                if (ix!=jx){
                    gmean += 2.0*xvalues(ix)*xvalues(jx)*Gijr(x,bp0,nmax,ix,jx);
                }
                else{
                    gmean += xvalues(ix)*xvalues(jx)*Gijr(x,bp0,nmax,ix,jx);
                } 
                //Get partials G(r)
                if(ix==0 && jx==0)          gtotal[1] = Gijr(x,bp0,nmax,ix,jx);
                else if(ix==0 && jx==(nsize+1)/2) gtotal[2] = Gijr(x,bp0,nmax,ix,jx);
                else if(ix==0 && jx==nsize-1)    gtotal[3] = Gijr(x,bp0,nmax,ix,jx);
                else if(ix==(nsize+1)/2 && jx==(nsize+1)/2) gtotal[4] = Gijr(x,bp0,nmax,ix,jx);
            }

        }
        gtotal[0] = gmean;
        return gtotal;
    }

    // Compute Gs real numbers
    template<class T>
    std::vector<T> Gs(T s, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        std::vector<T> gtotal = std::vector<T>(5,0.0);

        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                Qt(i,j) = Lmean*omega<T>(s+bp0,i,j);
                Qtm(i,j) = -Qt(i,j);
                if(i==j) Qtm(i,j) = 1.0+Qtm(i,j);

            }
        }

        //Invert matrix and all that
        Qt = Qt*(Qtm.inverse());
        T gmean = 0;
        double gmeanRE = 0;
        double gmeanIM = 0;

        double density = Density(bp);
        T dummy = 0.0;
        #pragma omp parallel for reduction(+:gmeanRE,gmeanIM) num_threads(16) private(dummy)
        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                dummy = sqrt(xvalues(i)*xvalues(j))*Qt(i,j);
                if (i!=j){
                    gmeanRE += std::real(2.0*dummy);
                    gmeanIM += std::imag(2.0*dummy);
                }
                else{
                    gmeanRE += std::real(dummy);
                    gmeanIM += std::imag(dummy);
                }
                //Get partials G(r)
                if(i==0 && j==0)          gtotal[1] = Qt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                else if(i==0 && j==(nsize+1)/2) gtotal[2] = Qt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                else if(i==0 && j==nsize-1)    gtotal[3] = Qt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                else if(i==(nsize+1)/2 && j==(nsize+1)/2) gtotal[4] = Qt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
            }
        }
        if constexpr (std::is_same_v<T, std::complex<double>>){
            gtotal[0] = std::complex<double>(gmeanRE/density, gmeanIM/density);
        }
        else{
            gtotal[0] = gmeanRE/density;
        }
        return gtotal;
    }

    // Compute Gijs, I need to change that, of course //Esto si todo va bien realmente no lo necesitare aunque lo puedo dejar por si acaso
    template<class T>
    T Gijs(T s,  double bp0, int ii, int jj){
        if (bp0 != bp)  SetPressure(bp0);
        T gtotal = 0;

        //Create matrix Qt
        Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> Qt;
        Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> Qtm;
        Qt.resize(nsize,nsize);
        Qtm.resize(nsize,nsize);
        #pragma omp parallel for num_threads(16)
        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                Qt(i,j) = Lmean*omega<T>(s+bp0,i,j);
                Qtm(i,j) = -Qt(i,j);
                if(i==j) Qtm(i,j) = 1.0+Qtm(i,j);

                // For symmetric indices too
                Qt(j,i) = Qt(i,j);
                Qtm(j,i) = Qtm(i,j);

            }
        }
        //Invert matrix and all that
        Qt = Qt*(Qtm.inverse());

        gtotal += Qt(ii,jj)/sqrt(xvalues(ii)*xvalues(jj));
        return gtotal/Density(bp);
    }

    //Compute Sk
    double Sk(double k, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double stotal = 0;

        return 1.0 + Density(bp0)*std::real(Gs<std::complex<double>>(std::complex<double>(0.0,k),bp0)[0]+Gs<std::complex<double>>(std::complex<double>(0.0,-k),bp0)[0]);

    }

    //Gr calculated from numerical Laplace inverse of G(s)
    std::vector<double> GrInv(double t, double bp0, int ntr, int meuler){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::vector<std::complex<double>> Gstemp = Gs<std::complex<double>>(xx,bp0);
        std::complex<double> suma = Gstemp[0]/2.0;
        std::complex<double> suma1 = Gstemp[1]/2.0;
        std::complex<double> suma2 = Gstemp[2]/2.0;
        std::complex<double> suma3 = Gstemp[3]/2.0;
        std::complex<double> suma4 = Gstemp[4]/2.0;

        std::complex<double> dummy;
        for (int i=1;i<ntr+1;i++){
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx, i*hh),bp0);
            dummy = std::complex<double>(pow(-1,i),0);
            suma += dummy * Gstemp[0];
            suma1 += dummy * Gstemp[1];
            suma2 += dummy * Gstemp[2];
            suma3 += dummy * Gstemp[3];
            suma4 += dummy * Gstemp[4];
        }

        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su1(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su2(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su3(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su4(meuler+2,std::complex<double>(0,0));
        su[0] = suma;
        su1[0] = suma1;
        su2[0] = suma2;
        su3[0] = suma3;
        su4[0] = suma4;

        //#avgsu
        //std::complex<double> avgsu(0,0);
        std::complex<double> avgsu(0,0);
        std::complex<double> avgsu1(0,0);
        std::complex<double> avgsu2(0,0);
        std::complex<double> avgsu3(0,0);
        std::complex<double> avgsu4(0,0);
        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            dummy = pow(-1,nn);
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx,nn*hh),bp0);
            su[i+1] = su[i] + dummy * Gstemp[0];
            su1[i+1] = su1[i] + dummy * Gstemp[1];
            su2[i+1] = su2[i] + dummy * Gstemp[2];
            su3[i+1] = su3[i] + dummy * Gstemp[3];
            su4[i+1] = su4[i] + dummy * Gstemp[4];
            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            //avgsu += bn*su[i];
            avgsu += bn*su[i+1];
            avgsu1 += bn*su1[i+1];
            avgsu2 += bn*su2[i+1];
            avgsu3 += bn*su3[i+1];
            avgsu4 += bn*su4[i+1];

        }
        double pp = pow(2.0,meuler);
        //std::complex<double> fun = uu*avgsu/(pp);
        std::vector<double> fun = std::vector<double>(5,0.0);
        fun[0] = real(uu*avgsu/(pp));
        fun[1] = real(uu*avgsu1/(pp));
        fun[2] = real(uu*avgsu2/(pp));
        fun[3] = real(uu*avgsu3/(pp));
        fun[4] = real(uu*avgsu4/(pp));

        return fun;
        //return 0.0;

    }

        //Gr calculated from numerical Laplace inverse of G(s)
    double GijrInv(double t, double bp0, int ntr, int meuler, int ii, int jj){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::complex<double> suma = Gijs(xx,bp0, ii, jj)/2.0;

        std::complex<double> argument(0,0);
        for (int i=1;i<ntr+1;i++){
            suma = suma + std::complex<double>(pow(-1,i),0) * Gijs(std::complex<double>(xx, i*hh),bp0, ii, jj);
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
            su[i+1] = su[i] + pow(-1,nn) * Gijs(std::complex<double>(xx,nn*hh),bp0, ii, jj);
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
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qt;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtm;

    template<class T>
    T omega(T s, int i, int j){
        return exp(-am(i,j)*s)/s;
    }

    template<class T>
    T omegap(T s, int i, int j){
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
        //for (int h=0; h<kvec.size()-1; h++){
        //    sa -= am(kvec.at(h),kvec.at(h+1));
        //}

        int h=0;
        while (sa >=0.0 && h<kvec.size()-1){
            sa -= am(kvec.at(h),kvec.at(h+1));
            h++;
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

    int nsize = 151; //Keep this value fixed.
    int nmax = 3; //Keep this value fixed.
   
    std::string tempvalue;
    //Ask for the excess pore size and make sure the value is valid.
    std::cout << "Excess pore size width [0-0.866] : ";
    std::getline(std::cin, tempvalue);
    double eps = std::stod(tempvalue);
    if(eps<=0.0 || eps>sqrt(3)/2){
        std::cout << "ERROR: Value of excess pore size parameter must be larger than 0 and lower than 0.866" << std::endl;
        return 1;
    }

    //Ask for the pressure value and make sure it is valid.
    std::cout << "Pressure (reduced units) : ";
    std::getline(std::cin, tempvalue);
    double bp = std::stod(tempvalue);
    if(bp<=0.0){
        std::cout << "ERROR: Value of pressure must be a positive number" << std::endl;
        return 1;
    }

    //Ask for what you want to compute.
    std::cout << "What do you want to compute? (1:Probability distribution, 2:Radial distribution function, 3:Structure factor) : ";
    std::getline(std::cin, tempvalue);
    int comValue = std::stoi(tempvalue);
    if(comValue!=1 && comValue!=2 && comValue!=3){
        std::cout << "ERROR: Wrong choice of number" << std::endl;
        return 1;
    }

    //Create an object
    Quasi1d f(nsize, eps);
    f.SetPressure(bp);


    double rmin, rmax, rcount;
    //Compute Sk
    if(comValue==3){
        
        rmin = 0;
        rmax = 12;
        rcount = 300;
        double rstep = (rmax-rmin)/rcount;

        double ri;
        std::vector<double> kvalues(rcount,0.0);
        std::vector<double> Skvalues(rcount,0.0);

  
        auto begin = std::chrono::high_resolution_clock::now();
        int bwidth = 70;
        double progress;
        double skmax = 0;
        double pkmax = 0;
        for(int k=0;k<=rcount;k++){
            ri = rmin + k*rstep;
            kvalues[k] = ri;
            Skvalues[k] = f.Sk(ri, bp);
            //std::cout << gvalues[k] << ", " << skmax << std::endl;
            if (Skvalues[k] > skmax){
                //std::cout << "Updating SkMax" << std::endl;
                skmax = Skvalues[k];
                pkmax = ri;
            }

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
             std::cout.flush();

        }
        std::cout << std::endl;

        std::cout << std::fixed << std::setprecision(5) << "The peak of the function is located at k = " << pkmax << std::endl;
    
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

        std::cout << "Elapsed time: " << elapsed.count()*1e-9 << std::endl;

        //Write results to file
        std::ofstream ofile;
        ofile.open("output_Sk.txt");
        for(int k=0; k<=rcount; k++){
            ofile << kvalues[k] << "\t" << Skvalues[k] << "\n";
        }
        ofile.close();
    }
    else if (comValue==2){    //Compute Gr
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
        std::vector<double> Grtemp;
        for(int k=0;k<=rcount;k++){
            ri = rmin + k*rstep;
            rvalues[k] = ri;
            //g(r)
            if (ri < (nmax+1)*f.amin){
                //std::cout << ri << " " << "Analytical" << std::endl;
                Grtemp = f.Gr(ri, bp, nmax);
                gvalues[k] = Grtemp[0];
                gvalues11[k] = Grtemp[1];
                gvalues13[k] = Grtemp[3];
                gvalues12[k] = Grtemp[2];
                gvalues22[k] = Grtemp[4];
            }
            else{
                //std::cout << ri << " " << "Numerical" << std::endl;
                Grtemp = f.GrInv(ri, bp, 100,30);
                gvalues[k] = Grtemp[0];
                gvalues11[k] = Grtemp[1];
                gvalues13[k] = Grtemp[3];
                gvalues12[k] = Grtemp[2];
                gvalues22[k] = Grtemp[4];           
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
            std::cout.flush();

        }
        std::cout << std::endl;
    
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

        std::cout << "Elapsed time: " << elapsed.count()*1e-9 << std::endl;

        //Write results to file
        std::ofstream ofile;
        ofile.open("output_Gr.txt");
        for(int k=0;k<=(int)rcount;k++){
            ofile << rvalues[k] << "\t" << gvalues[k] << " " << gvalues11[k] << " " << gvalues12[k]<< " " << gvalues13[k]<< " " << gvalues22[k] << "\n";
            //ofile << rvalues[k] << "\t" << gvalues[k] << "\n";
        }
        ofile.close();
    }

  return 0;
}