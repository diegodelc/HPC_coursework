#include <iostream>
using namespace std;
#include <cmath>
#include <cblas.h>

#include <iomanip>
#include <fstream>

#include "ShallowWater.h"


int main(int argc, char **argv)
{
    // define parameters
    /*
    --dt arg    //Time-step to use.
    --T arg     //Total integration time.
    --Nx arg    //Number of grid points in x
    --Ny arg    //Number of grid points in y
    --ic arg    //Index of the initial condition to use (1-4)
    */
    double dt = 0.1;
    double T = 20;
    int Nx = 300;
    int Ny = 300;
    int ic = 4;
    
    //Instantiate class and parameters via constructor
    int whichIntegrationMethod = 1; // 1: For loop implementation,  2: BLAS implementation
    
    ShallowWater myInstance(dt,T,Nx,Ny,ic,whichIntegrationMethod);
    
    // Initialise simulation
    myInstance.SetInitialConditions();
    
    // Perform time integration
    myInstance.TimeIntegrate();
    
    
    // output to file
    cout << endl << "Writing output to file: ";
    
    ofstream vOut("output.txt", ios::out | ios::trunc);
    vOut.precision(5);
    for (int yInd = 0;yInd<Ny;yInd++) {
        for (int xInd = 0;xInd<Nx;xInd++) {
        vOut << xInd*myInstance.dx << " "
            << yInd*myInstance.dx << " "
            << myInstance.yn[yInd*Ny+xInd] << " "
            << myInstance.yn[yInd*Ny+xInd + Nx*Ny] << " "
            << myInstance.yn[yInd*Ny+xInd + 2*Nx*Ny] << " " << endl;
        }
    }
    vOut.close();
    cout << "DONE" << endl;
}
