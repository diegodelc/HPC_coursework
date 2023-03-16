#include <iostream>
using namespace std;
#include <cmath>

class ShallowWater {
public:
        double dt;
        double T;
        int Nx;
        int Ny;
        int ic;
        
        double* u;
        double* v;
        double* h;
        
        double dx = 1;
        
        ShallowWater(double dt_in,double T_in,
                            int Nx_in,int Ny_in,
                            int ic_in) {
            cout << "SetParameters" << endl;
            dt = dt_in;
            T = T_in;
            Nx = Nx_in;
            Ny = Ny_in;
            ic = ic_in;
            
            u = new double[Nx*Ny];
            v = new double[Nx*Ny];
            h = new double[Nx*Ny];
        };
        
        void SetInitialConditions() {
            /*
            This function initialises the u, v and h arrays to the values defined in the assignement brief
            
                u and v are always zero, so they are initialised separately from h
                h is initialised as one of four intial conditions, specified via the parameter 'ic'
            */
            
            cout << "SetInitialConditions" << endl;
            
            
            for (int x = 0;x<Nx; x++) {
                for (int y = 0;y<Ny; y++) {
                    u[x*Ny + y] = 0;
                    v[x*Ny + y] = 0;
                }
            }
            

            //  This is slower than it has to be, it is evaluating the if Nx*Ny times, 
            //  would be better to evaluate if and then go into for
            for (int x = 0;x<Nx; x++) {
                for (int y = 0;y<Ny; y++) {
                    if (ic == 1) {
                        h[x*Ny + y] = initialCond1(x*dx,y*dx);
                    } else if (ic == 2) {
                        h[x*Ny + y] = initialCond2(x*dx,y*dx);
                    } else if (ic == 3) {
                        h[x*Ny + y] = initialCond3(x*dx,y*dx);
                    } else if (ic == 4) {
                        h[x*Ny + y] = initialCond4(x*dx,y*dx);
                    }
                }
            }
            
            
        };
        void TimeIntegrate() {
            cout << "TimeIntegrate" << endl;
            // time propagation (for or while)
    
                // calculate k1, k2, k3, k4
                
                //find yn+1
                
                //t += dt
            
        };
private:
    double initialCond1(double x,double y) {
        return 10 + exp(-(x-50)*(x-50)/25);
    };
    double initialCond2(double x,double y) {
        return 10 + exp(-(y-50)*(y-50)/25);
    };
    double initialCond3(double x,double y) {
        return 10 + exp(-((x-50)*(x-50) + (y-50)*(y-50))/25);
    };
    double initialCond4(double x,double y) {
        return 10 + exp(-(x-25)*(x-25)/25) + exp(-(x-75)*(x-75)/25);
    };
    
};


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
    double T = 5;
    int Nx = 50;
    int Ny = 30;
    int ic = 1;
    
    //Instantiate class and parameters via constructor
    ShallowWater myInstance(dt,T,Nx,Ny,ic);
    cout << myInstance.dt << endl;
    cout << myInstance.T << endl;
    cout << myInstance.Nx << endl;
    cout << myInstance.Ny << endl;
    cout << myInstance.ic << endl;
    
    
    // initial conditions and grids
    myInstance.SetInitialConditions();
    
    // perform time integration
    myInstance.TimeIntegrate();
    
    
    
    
    
    // output to file
}
