#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif 

#ifndef QUASI1D_H_INCLUDED
#define QUASI1D_H_INCLUDED

//Efficient way of computing binomial coefficients
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

//Class Quasi1d
class Quasi1d {       
  public:             
    int nsize; //size of the discretization (number of species).
    double eps; //excess width of the pore
    double bp; //pressure
    double amin; //min distance between two particles

    double Lmean; //Lmean coming from the eigenvalue equation.
    Eigen::Matrix<double, Eigen::Dynamic, 1> xvalues; //vector of molar fractions.

    int nprocs = omp_get_num_procs();

    //Constructor
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

    //Set pressure of the system and update xvalues and Lmean.
    int SetPressure(double bp0){
        //Create matrix M
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                M(i,j) = omega(bp0,i,j);
            }
        }

        //Solve eigenvalue equation
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
        this->bp = bp0;
        return 0;

    }

    //If instead of pressure the user wishes to set the density, use this function.
    void SetDensity(double rho){
        //Bolzano method to compute the pressure associated with that density
        double a = 0.000001;
        double b=170;
  
        double c = a;
        while ((b-a) >= 0.0000001){
 
            // Find middle point
            c = (a+b)/2.0;
           
            // Decide the side to repeat the steps
            if ((Density(c)-rho)*(Density(a)-rho) < 0)   b = c;
            else  a = c;
        }

        if (abs(c-170.0)<=0.0001){
            std::cout << "WARNING: required density set too high, numerical errors might occur. Setting new density to " << Density(c) << std::endl;
        }

        //Set the corresponding pressure
        SetPressure(c);


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

    // Compute Gij(x) analiticamente para un nmax
    double Gijx(double x, double bp0, int nmax, int i, int j){
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

    // Compute  G(x)
    std::vector<double> Gx(double x, double bp0, int nmax){
        if (bp0 != bp)  SetPressure(bp0);
        std::vector<double> gtotal = std::vector<double>(5,0.0);

        //Optimizacion: vemos en que capa de proximos vecinos cae x y ajustamos el valor de nmax si es necesario
        int nshell = int(x/amin);
        if (nshell == 0) return gtotal;
        else if (nshell < nmax)  nmax = nshell; //keep lowest n value because they are going to yield the same result

        double gtemp = 0.0;
        double gmean = 0.0;
        
        //#pragma omp parallel for reduction(+:gmean) num_threads(16) schedule(dynamic)
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                //gtemp = Gijx(x,bp0,nmax,ix,jx);
                if (ix!=jx){
                    gmean += 2.0*xvalues(ix)*xvalues(jx)*Gijx(x,bp0,nmax,ix,jx);
                }
                else{
                    gmean += xvalues(ix)*xvalues(jx)*Gijx(x,bp0,nmax,ix,jx);
                } 
                //Get partials G(r)
                if(ix==0 && jx==0)          gtotal[1] = Gijx(x,bp0,nmax,ix,jx);
                else if(ix==0 && jx==(nsize+1)/2) gtotal[2] = Gijx(x,bp0,nmax,ix,jx);
                else if(ix==0 && jx==nsize-1)    gtotal[3] = Gijx(x,bp0,nmax,ix,jx);
                else if(ix==(nsize+1)/2 && jx==(nsize+1)/2) gtotal[4] = Gijx(x,bp0,nmax,ix,jx);
            }

        }
        gtotal[0] = gmean;
        return gtotal;
    }

    // Compute G(s)
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
        #pragma omp parallel for reduction(+:gmeanRE,gmeanIM) num_threads(nprocs) private(dummy)
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

    // Compute Gij(s)
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

    //Compute S(k)
    double Sk(double k, double bp0, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        if (bp0 != bp)  SetPressure(bp0);
        double stotal = 0;

        return 1.0 + Density(bp0)*std::real(Gs<std::complex<double>>(std::complex<double>(0.0,k),bp0,Qtt, Qttm)[0]+Gs<std::complex<double>>(std::complex<double>(0.0,-k),bp0, Qtt, Qttm)[0]);

    }

    //G(x) calculated from numerical Laplace inverse of G(s)
    std::vector<double> GxInv(double t, double bp0, int ntr, int meuler, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::vector<std::complex<double>> Gstemp = Gs<std::complex<double>>(xx,bp0,Qtt, Qttm);
        std::complex<double> suma = Gstemp[0]/2.0;
        std::complex<double> suma1 = Gstemp[1]/2.0;
        std::complex<double> suma2 = Gstemp[2]/2.0;
        std::complex<double> suma3 = Gstemp[3]/2.0;
        std::complex<double> suma4 = Gstemp[4]/2.0;

        std::complex<double> dummy;
        for (int i=1;i<ntr+1;i++){
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx, i*hh),bp0,Qtt, Qttm);
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
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx,nn*hh),bp0,Qtt, Qttm);
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

    //Gij(x) calculated from numerical Laplace inverse of G(s)
    double GijxInv(double t, double bp0, int ntr, int meuler, int ii, int jj){
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

        #pragma omp parallel for reduction(+:gtotal) num_threads(nprocs) schedule(dynamic)
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
    
    //Compute G(r)
    double Gr(double x, double bp0, int nmax, bool b2D=0){
        if (bp0 != bp)  SetPressure(bp0);
        double gtotal = 0.0;

        //Optimizacion: vemos en que capa de proximos vecinos cae x y ajustamos el valor de nmax si es necesario
        int nshell = int(x/amin);
        if (nshell == 0) return gtotal;
        else if (nshell < nmax)  nmax = nshell; //keep lowest n value because they are going to yield the same result
        
        //#pragma omp parallel for reduction(+:gtotal) num_threads(16) schedule(dynamic)
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                if (ix!=jx){
                    gtotal += 2.0*xvalues(ix)*xvalues(jx)*Gijx(sqrt(x*x-(y(ix)-y(jx))*(y(ix)-y(jx))),bp0,nmax,ix,jx);
                }
                else{
                    gtotal += xvalues(ix)*xvalues(jx)*Gijx(sqrt(x*x-(y(ix)-y(jx))*(y(ix)-y(jx))),bp0,nmax,ix,jx);
                } 
            }

        }
        return gtotal;
    }
   
   //G(r) calculated from numerical Laplace inverse of G(s)
   double GrInv(double x, double bp0, int ntr, int meuler){
        if (bp0 != bp)  SetPressure(bp0);
        double gtotal = 0.0;

        double dy2 = (eps/(nsize-1))*(eps/(nsize-1));

        //#pragma omp parallel for reduction(+:gtotal) num_threads(16) schedule(dynamic)
        for (int k=0; k<nsize; k++){
                Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtt;
                Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qttm;
                Qtt.resize(nsize,nsize);
                Qttm.resize(nsize,nsize);
                gtotal += GrkiInv(sqrt(x*x-dy2*k*k), bp0, k, ntr, meuler,Qtt,Qttm);
                //std::cout << k << ", " << sqrt(x*x-dy2*k*k)<< ", " << bp0 << ","<<ntr << "," << meuler << ","  << gtotal << std::endl;
        }
        return gtotal;
    }

    //Gk(r) calculated from numerical Laplace inverse of G(s)
    double GrkiInv(double t, double bp0, int k, int ntr, int meuler, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::complex<double> Gstemp = Gkis<std::complex<double>>(xx,bp0,k,Qtt,Qttm);
        std::complex<double> suma = Gstemp/2.0;

        std::complex<double> dummy;
        for (int i=1;i<ntr+1;i++){
            Gstemp = Gkis<std::complex<double>>(std::complex<double>(xx, i*hh),bp0,k, Qtt, Qttm);
            dummy = std::complex<double>(pow(-1,i),0);
            suma += dummy * Gstemp;
        }
        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        su[0] = suma;

        //#avgsu
        //std::complex<double> avgsu(0,0);
        std::complex<double> avgsu(0,0);

        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            dummy = pow(-1,nn);
            Gstemp = Gkis<std::complex<double>>(std::complex<double>(xx,nn*hh),bp0,k, Qtt, Qttm);
            su[i+1] = su[i] + dummy * Gstemp;

            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            //avgsu += bn*su[i];
            avgsu += bn*su[i+1];

        }
        double pp = pow(2.0,meuler);
        //std::complex<double> fun = uu*avgsu/(pp);
        double fun = 0.0;
        fun = real(uu*avgsu/(pp));
        return fun;
        //return 0.0;

    }

    template<class T>
    T Gkis(T s, double bp0, int kidx, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        if (bp0 != bp)  SetPressure(bp0);
        T gtotal = 0.0;

        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                //std::cout << i << " " << j << std::endl;
                Qtt(i,j) = Lmean*omega<T>(s+bp0,i,j);
                Qttm(i,j) = -Qtt(i,j);
                if(i==j) Qttm(i,j) = 1.0+Qttm(i,j);
            }
        }

        //Invert matrix and all that
        Qtt = Qtt*(Qttm.inverse());
        double density = Density(bp);
        for (int i=0; i<nsize-kidx; i++){
                gtotal += sqrt(xvalues(i)*xvalues(i+kidx))*Qtt(i,i+kidx)/density; //multiplicado por xi*xi y dividido entre sqrt(xi xj) da la simplificacion que sale aqui
        }
        if (kidx>0) gtotal = 2.0*gtotal;
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

        if (sa >= 0){
            return exp(-bp0*x)*pow(sa, n-1.0)/tgamma(n); //tgamma(n+1)=factorial(n)
        }
        else{
            return 0.0;
        }
    }
};

#endif
