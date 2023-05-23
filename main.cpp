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

#include "quasi1d.hpp"
#include "iofunctions.hpp"


int main(int argc, char* argv[]) {

    InputData indata;
    if(ReadInputFile("input.dat", &indata)){
        std::cout << "Exiting...\n" << std::endl;
        std::cout << "\n Press enter to exit \n" << std::endl;
        std::cin.get();
        return 1;
    }
    if(CheckData(&indata)){
        std::cout << "Exiting...\n" << std::endl;
        std::cout << "\n Press enter to exit \n" << std::endl;
        std::cin.get();
        return 1;
    }

    int nsize = 251; //Keep this value fixed.
    int nmax = 3; //Keep this value fixed.
   
    std::string tempvalue;
    //Create an object
    Quasi1d f = Quasi1d(nsize, indata.eps);
    if (indata.density>0.0){
        f.SetDensity(indata.density);
        std::cout << "The pressure of the system has been set to " << f.bp << std::endl;
        indata.bp = f.bp;
    }
    else{
        f.SetPressure(indata.bp);
        std::cout << "The density of the system has been set to " << f.Density(indata.bp) << std::endl;
    }

    //vectors to store results
    double xi;
    double rstep = (indata.xmax-indata.xmin)/(indata.npoints-1);
    std::vector<double> rvalues(indata.npoints,0.0);
    std::vector<double> Fvalues(indata.npoints,0.0);

    //Compute Probability
    if(indata.comValue==1){
        //Ask for what you want to compute.
        std::cout << "Which n-neighbor distribution do you want to compute (1, 2 or 3)?: ";
        std::getline(std::cin, tempvalue);
        int nn = std::stoi(tempvalue);
        if(nn!=1 && nn!=2 && nn!=3){
            std::cout << "ERROR: Please choose 1, 2 or 3" << std::endl;
            return 1;
        }

        double normalization = 0.0; //To check normalization
        std::cout << "Computing..." <<std::endl;
        for(int k=0;k<indata.npoints;k++){
            xi = indata.xmin + k*rstep;
            rvalues[k] = xi;
            //Probability Density
            Fvalues[k] = f.ProbDensity(xi,nn,f.bp);
            normalization += rstep*Fvalues[k];
             //Progress bar
            ProgressBar(double(k)/(double(indata.npoints)-1),70);

        }
        std::cout << std::endl;

        std::cout << "Check Normalization: " << normalization << std::endl;
        std::cout << "Done!\n" << std::endl;
        //Write results to file
        std::cout << "******************\nWriting results to file " << "output_Pn.txt" << std::endl;
        std::cout << "Formatting:\n x \t Pn(x)\n\n";
        WriteToFile("output_Pn.txt", indata.npoints,&rvalues,&Fvalues);
        std::cout << "******************" << std::endl;

    }
    else if (indata.comValue==2){    //Compute Gr
        std::vector<double> Grtemp;
        std::vector<double> gvalues11(indata.npoints,0.0);
        std::vector<double> gvalues12(indata.npoints,0.0);
        std::vector<double> gvalues22(indata.npoints,0.0);
        std::vector<double> gvalues13(indata.npoints,0.0);
    
        auto begin = std::chrono::high_resolution_clock::now();
        
        std::cout << "Computing..." << std::endl;
        for(int k=0;k<indata.npoints;k++){
            xi = indata.xmin + k*rstep;
            rvalues[k] = xi;
            //g(r)
            if (xi < (nmax+1)*f.amin){
                Grtemp = f.Gx(xi, f.bp, nmax);
                Fvalues[k] = Grtemp[0];
                gvalues11[k] = Grtemp[1];
                gvalues13[k] = Grtemp[3];
                gvalues12[k] = Grtemp[2];
                gvalues22[k] = Grtemp[4];
            }
            else{
                Grtemp = f.GxInv(xi, f.bp, 100,30);
                Fvalues[k] = Grtemp[0];
                gvalues11[k] = Grtemp[1];
                gvalues13[k] = Grtemp[3];
                gvalues12[k] = Grtemp[2];
                gvalues22[k] = Grtemp[4];           
            }

            //Progress bar
            ProgressBar(double(k)/(double(indata.npoints)-1),70);

        }
        std::cout << std::endl;
        std::cout << "Done!\n" << std::endl;

        //Write results to file
        std::cout << "******************\nWriting results to file " << "output_Gx.txt" << std::endl;
        std::cout << "Formatting:\n x \t G(x)\n\n";
        std::ofstream ofile;
        ofile.open("output_Gx.txt");
        for(int k=0;k<(int)indata.npoints;k++){
            ofile << rvalues[k] << "\t" << Fvalues[k] << " " << gvalues11[k] << " " << gvalues12[k]<< " " << gvalues13[k]<< " " << gvalues22[k] << "\n";
        }
        ofile.close();
        std::cout << "******************" << std::endl;
    }
    else if(indata.comValue==3){
         
        double skmax = 0;
        double pkmax = 0;
        std::cout << "Computing..." <<std::endl;
        for(int k=0; k<indata.npoints; k++){
            xi = indata.xmin + k*rstep;
            rvalues[k] = xi;
            Fvalues[k] = f.Sk(xi, f.bp);

            if (Fvalues[k] > skmax){
                skmax = Fvalues[k];
                pkmax = xi;
            }

             //Progress bar
            ProgressBar(double(k)/(double(indata.npoints)-1),70);

        }
        std::cout << std::endl;

        std::cout << "Done!\n" << std::endl;

        //std::cout << std::fixed << std::setprecision(9) << "Peak of S()q= " << f.Density(f.bp) << " " << pkmax << std::endl;

        //Write results to file
        std::cout << "******************\nWriting results to file " << "output_Sq.txt" << std::endl;
        std::cout << "Formatting:\n q \t S(q)\n\n";
        WriteToFile("./output_Sq.txt", indata.npoints, &rvalues, &Fvalues);
        std::cout << "******************" << std::endl;

    }

    //std::cout << "\n Press enter to exit \n" << std::endl;
    //std::cin.get();

  return 0;
}