#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

struct InputData{
    double eps = 0.0;
    double bp = 0.0;
    double density = 0.0;
    int comValue = 1;
    double xmin = 0.0;
    double xmax = 1.0;
    int npoints = 100;
};

//Function definitions
int ReadInputFile(std::string, InputData*);
int CheckData(InputData*);
void WriteInputFileText();
void WriteToFile(std::string, int, std::vector<double>*, std::vector<double>*);
void ProgressBar(double, int);



int ReadInputFile(std::string infile, InputData *indata){
    //Read input file
    std::ifstream inFile;
    //Open file
    inFile.open(infile);
    if (!inFile.is_open()){
        WriteInputFileText();
        std::cout << "Could not file input.dat. A new standard input file has been generated to work with. Please use that file to change parameters of the system" << std::endl;
        return 1;
    }
    int lineNum = 0;
    std::string line;
    while (std::getline(inFile, line))
    {
        std::istringstream iss(line);
        if (line.length() > 0 && line[0] == '#') continue;
        std::string v1;
        double v2;
        iss >> v1 >> v2;
        try{
            if (v1=="eps:")       indata->eps = double(v2);
            if (v1=="pressure:")        indata->bp = double(v2);
            if (v1=="density:")        indata->density = v2;
            if (v1=="quantity:")        indata->comValue = int(v2);
            if (v1=="min:")        indata->xmin = v2;
            if (v1=="max:")        indata->xmax = v2;
            if (v1=="npoints:")        indata->npoints = int(v2);
        }
        catch(std::exception &e){
            std::cout << "Something went wrong when reading input file. Please make sure everything is correct." << std::endl;
            WriteInputFileText();
            return 1;
        }
    }
    return 0;
}

int CheckData(InputData *indata){
    std::cout << "*************************\n SYSTEM INFO:" <<std::endl;
    double density_max = 1.0/sqrt(1.0-indata->eps*indata->eps);

    //Check value of epsilon
    if(indata->eps<=0.0 || indata->eps>sqrt(3)/2){
        std::cout << "eps = " << indata->eps << std::endl;
        std::cout << "ERROR: Value of excess pore size parameter must be larger than 0 and lower than 0.866" << std::endl;
        return 1;
        }
    else{
        std::cout << "eps:" << indata->eps << std::endl;
    }

    if(indata->density<0.0){
        std::cout << "ERROR: Value of density must be positive (overwrites pressure value) or zero (it is not read)" << std::endl;
        return 1;
    }
    else if (indata->density>density_max) {
        std::cout << "ERROR: Value of density must be lower than " << density_max << std::endl;
        return 1;
    }
    else if (indata->density>0.0) {
        std::cout << "density:" << indata->density << " (overwrites pressure value)" << std::endl;
    }
    else if (indata->density == 0.0){
        if (indata->bp < 0.0){
            std::cout << "ERROR: Value of pressure must be a positive number" << std::endl;
            return 1;
        }
        else{
            std::cout << "pressure:" << indata->bp << std::endl;
    
        }
    }
    if(indata->comValue!=1 && indata->comValue!=2 && indata->comValue!=3){
        std::cout << "ERROR: Wrong choice of number" << std::endl;
        return 1;
        }
        else{
            std::cout << "Compute quantity:" << indata->comValue << std::endl;
    }
    if(indata->xmin<0){
        std::cout << "ERROR: Value of xmin must be a positive number." << std::endl;
        return 1;
        }
        else{
            std::cout << "min:" << indata->xmin << std::endl;
    }
    if(indata->xmax < indata->xmin){
            std::cout << "ERROR: max must be bigger than min." << std::endl;
            return 1;
        }
        else{
            std::cout << "max:" << indata->xmax << std::endl;
    }
    if(indata->npoints<=1){
        std::cout << "ERROR: number of points must be 2 or higher." << std::endl;
        return 1;
        }
        else{
            std::cout << "npoints:" << indata->npoints << std::endl;
    }
    std::cout << "*************************\n" <<std::endl;
    return 0;
}

void WriteInputFileText(){
    std::ofstream ofile;
    ofile.open("./input.dat");
    ofile << "#Excess pore size width [0-0.866]\n";
    ofile << "eps: 0.866\n\n";

    ofile << "#Set pressure or density (setting density>0 will overwrite pressure value)\n";
    ofile << "pressure: 1.0 \n";
    ofile << "density: 0.0\n\n";

    ofile << "#Quantity to compute (1:Probability distribution, 2:Radial distribution function, 3:Structure factor)\n";
    ofile << "quantity: 2\n\n";

    ofile << "#Range in the form min max npoints\n";
    ofile << "min: 0\n";
    ofile << "max: 5\n";
    ofile << "npoints: 100\n";
    ofile.close();
}

void WriteToFile(std::string outfile, int npoints, std::vector<double> *xvalues, std::vector<double> *Fvalues){
            //Write results to file
        std::ofstream ofile;
        ofile.open(outfile);
        for(int k=0; k<npoints; k++){
            ofile << xvalues->at(k) << "\t" << Fvalues->at(k) << "\n";
        }
        ofile.close();
}

void ProgressBar(double progress, int bwidth){
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

  #endif