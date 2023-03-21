#include <iostream>
using namespace std;
#include <cmath>
#include <cblas.h>

#include <iomanip>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <omp.h>

#include "ShallowWater.h"

namespace po = boost::program_options;

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
    //double dt = 0.1;
    //double T = 20;
    //int Nx = 100;
    //int Ny = 100;
    //int ic = 4;
    
        // Boost program options
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "Prints list of options")
        ("dt", po::value<double>()->default_value(0.1), "Time step")
        ("T", po::value<double>()->default_value(5.0), "Total simulation time")
        ("Nx", po::value<int>()->default_value(100), "Number of grid points in x")
        ("Ny", po::value<int>()->default_value(100), "Number of grid points in y")
        ("ic", po::value<int>()->default_value(4), "Initial conditions");
        

    // Default execution:
    //     ./main --dt 0.1 --T 5 --Nx 100 --Ny 100 --ic 4

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << "\n";
        return 1;
    }

    // Loading parameters into memory
    double dt     = vm["dt"].as<double>();
    double T     = vm["T"].as<double>();
    int Nx     = vm["Nx"].as<int>();
    int Ny     = vm["Ny"].as<int>();
    int ic     = vm["ic"].as<int>();
    
    
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
    
    return 0;
}
