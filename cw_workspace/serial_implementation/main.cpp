#include <iostream>
using namespace std;
#include <cmath>
#include <cblas.h>

#include <iomanip>
#include <fstream>

class ShallowWater {
public:
        double dt;
        double T;
        int Nx;
        int Ny;
        int ic;
        
        double* yn;
        
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
            yn = new double[3*Nx*Ny];
            
            //u = new double[Nx*Ny];
            //v = new double[Nx*Ny];
            //h = new double[Nx*Ny];
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
                    yn[x*Ny + y] = 0; //u
                    yn[Nx*Ny + x*Ny + y] = 0; //v
                }
            }
            

            //  This is slower than it has to be, it is evaluating the if Nx*Ny times, 
            //  would be better to evaluate if and then go into for
            for (int x = 0;x<Nx; x++) {
                for (int y = 0;y<Ny; y++) {
                    if (ic == 1) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond1((double)x*dx,(double)y*dx);
                    } else if (ic == 2) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond2((double)x*dx,(double)y*dx);
                    } else if (ic == 3) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond3((double)x*dx,(double)y*dx);
                    } else if (ic == 4) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond4((double)x*dx,(double)y*dx);
                    }
                }
            }
            
            
            
        };
        
        void TimeIntegrate() {
            cout << endl << "TimeIntegrate" << endl;
            double* dudx = new double[Nx*Ny];
            double* dudy = new double[Nx*Ny];
            
            double* dvdx = new double[Nx*Ny];
            double* dvdy = new double[Nx*Ny];
            
            double* dhdx = new double[Nx*Ny];
            double* dhdy = new double[Nx*Ny];
            
            
            
            double* temp = new double[3*Nx*Ny];
            double* allKs = new double[3*Nx*Ny];
            // time propagation (for or while)
            int counter = 0;
            for (double time = 0; time < T; time += dt) { //T and dt are double, so easier to make time double than try and cast them to int
                
                
                
                
                //k1 = calcF(yn)
                calcF(  yn,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);                
                cblas_dcopy(3*Nx*Ny,temp,1,allKs,1);
                
                //k2 = calcF(yn + dt*k1/2);
                multVectByConst(3*Nx*Ny,temp,dt/2,temp,'R');
                addTwoVectors(3*Nx*Ny,yn,temp,1,temp);
                calcF(  temp,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);
                addTwoVectors(3*Nx*Ny,allKs,temp,2,allKs);
                
                //k3 = calcF(yn + dt*k2/2);
                multVectByConst(3*Nx*Ny,temp,dt/2,temp,'R');
                addTwoVectors(3*Nx*Ny,yn,temp,1,temp);
                calcF(  temp,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);
                addTwoVectors(3*Nx*Ny,allKs,temp,2,allKs);
                
                //k4 = calcF(yn + dt*k3);
                multVectByConst(3*Nx*Ny,temp,dt,temp,'R');
                addTwoVectors(3*Nx*Ny,yn,temp,1,temp);
                calcF(  temp,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);
                addTwoVectors(3*Nx*Ny,allKs,temp,1,allKs);
                
                
                //yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt;
                multVectByConst(3*Nx*Ny,allKs,dt/6,allKs,'R');
                addTwoVectors(3*Nx*Ny,yn,allKs,1,yn);//inefficient, mult by zero inneccessary
                
                
                
                cout << "Finished iteration " << time/dt << "/" << T/dt << endl;
                if (counter == 1) {
                    for (int i=0;i<Nx*Ny;i++) {
                        cout << dvdx[2*Nx*Ny + i] << endl;
                    }
                }
                
                counter++;
                
                
                
            };

            delete[] dudx;
            delete[] dudy;
            
            delete[] dvdx;
            delete[] dvdy;
            
            delete[] dhdx;
            delete[] dhdy;
            
            delete[] allKs;
            delete[] temp;
        };
        
private:
    double stencil[7] = {-0.0167, 0.1500, -0.7500, 0, 0.7500, -0.1500, 0.0167};
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
    
    void multVectByConst(const int& n,const double* vect1,const double& constVal,double* ans,const char& HOW) {
        /*
         HOW: 'A' : add to answer
              'R' : replace answer (destroy previous values)
          */
        if (HOW == 'R') {
            for (int i=0; i<n; i++) {
                ans[i] = constVal*vect1[i];
            }
        } else if (HOW == 'A') {
            for (int i=0; i<n; i++) {
                ans[i] += constVal*vect1[i];
            }
        } else { //raiseAssertion
            cout << "Error: Option < " << HOW << "> not Implemented, only add (A) or replace (R)" << endl;
        }
    };
    
    void addTwoVectors(const int& n,const double* vect1,const double* vect2,const double& sign,double* ans) {
        for (int i = 0; i<n; i++) {
            ans[i] = vect1[i] + sign * vect2[i];
        }
    };
    
    void multVectByVect(const int& n,const double* vect1,const double* vect2,const double& sign,double* ans,const char& HOW) {
        /*
         HOW: 'A' : add to answer
              'R' : replace answer (destroy previous values)
          */
        if (HOW == 'R') {
            for (int i=0; i<n; i++) {
                ans[i] = sign*vect1[i]*vect1[i];
            }
        } else if (HOW == 'A') {
            for (int i=0; i<n; i++) {
                ans[i] += sign*vect1[i]*vect1[i];
            }
        } else { //raiseAssertion
            cout << "ERROR: Option < " << HOW << "> not Implemented, only add (A) or replace (R)" << endl;
        }
    };
    
    void calcF( double* yn,
                double* dudx,double* dudy,
                double* dvdx,double* dvdy,
                double* dhdx,double* dhdy,
                double* f) {
        
        //cblas_dcopy(Nx*Ny,yn,1,u,1);
        //cblas_dcopy(Nx*Ny,yn+Nx*Ny,1,v,1);
        //cblas_dcopy(Nx*Ny,yn+2*Nx*Ny,1,h,1);
        
        
        //NEED TO CHANGE NAMES TO U_POS SO THEY DON'T OVERWRITE THE ATTRIBUTES
        //ENSURE WE ARE NOT OVERWRITING YN, WE NEED THIS INTACT FOR EACH K EVALUATION
        double* u_pos = yn;
        double* v_pos = yn+Nx*Ny;
        double* h_pos = yn+2*Nx*Ny;
        
        
        derXFor(u_pos,dudx);
        derYFor(u_pos,dudy);
        
        derXFor(v_pos,dvdx);
        derYFor(v_pos,dvdy);
        
        derXFor(h_pos,dhdx);
        derYFor(h_pos,dhdy);
        
        //f1 = - (g*dhdx + u.*dudx) - (v.*dudy);
        multVectByConst(Nx*Ny,dhdx,-9.81,f,'R');
        multVectByVect(Nx*Ny,dudx,u_pos,-1,f,'A');
        multVectByVect(Nx*Ny,dudy,v_pos,-1,f,'A');
        
        //f2 = - (u.*dvdx) - (g*dhdy + v.*dvdy);
        multVectByVect(Nx*Ny,u_pos,dvdx,-1,f+Nx*Ny,'R');
        multVectByConst(Nx*Ny,dhdy,-9.81,f+Nx*Ny,'A');
        multVectByVect(Nx*Ny,v_pos,dvdy,-1,f+Nx*Ny,'A');
        
        //f3 = - (u.*dhdx + h.*dudx) - (v.*dhdy + h.*dvdy);
        multVectByVect(Nx*Ny,u_pos,dhdx,-1,f+2*Nx*Ny,'R');
        multVectByVect(Nx*Ny,h_pos,dudx,-1,f+2*Nx*Ny,'A');
        multVectByVect(Nx*Ny,v_pos,dhdy,-1,f+2*Nx*Ny,'A');
        multVectByVect(Nx*Ny,h_pos,dvdy,-1,f+2*Nx*Ny,'A');
    };
    
    void derXFor(double* data, double* derivative) {
        double* temp = new double[7];
        for (int col = 0; col < Nx; col++) {
            for (int row = 0; row < Ny; row++) {
                if (3<=col && col<=Nx-4) {
                    temp[0] = data[row*Ny+col-3];
                    temp[1] = data[row*Ny+col-2];
                    temp[2] = data[row*Ny+col-1];
                    temp[3] = data[row*Ny+col]; //multiplied by zero
                    temp[4] = data[row*Ny+col+1];
                    temp[5] = data[row*Ny+col+2];
                    temp[6] = data[row*Ny+col+3];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (col == 0) { //first column
                    temp[0] = data[row*Ny+Ny-3];
                    temp[1] = data[row*Ny+Ny-2];
                    temp[2] = data[row*Ny+Ny-1];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[row*Ny+col+1];
                    temp[5] = data[row*Ny+col+2];
                    temp[6] = data[row*Ny+col+3];

                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (col == 1) { //second column
                    temp[0] = data[row*Ny+Ny-2];
                    temp[1] = data[row*Ny+Ny-1];
                    temp[2] = data[row*Ny+col-1];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[row*Ny+col+1];
                    temp[5] = data[row*Ny+col+2];
                    temp[6] = data[row*Ny+col+3];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (col == 2) { //third column
                    temp[0] = data[row*Ny+Ny-1];
                    temp[1] = data[row*Ny+col-2];
                    temp[2] = data[row*Ny+col-1];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[row*Ny+col+1];
                    temp[5] = data[row*Ny+col+2];
                    temp[6] = data[row*Ny+col+3];

                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (col == Nx-1) { //last column
                    temp[0] = data[row*Ny+col-3];
                    temp[1] = data[row*Ny+col-2];
                    temp[2] = data[row*Ny+col-1];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[row*Ny+0];
                    temp[5] = data[row*Ny+1];
                    temp[6] = data[row*Ny+2];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (col == Nx-2) { //second-to-last column
                    temp[0] = data[row*Ny+col-3];
                    temp[1] = data[row*Ny+col-2];
                    temp[2] = data[row*Ny+col-1];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[row*Ny+col+1];
                    temp[5] = data[row*Ny+0];
                    temp[6] = data[row*Ny+1];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (col == Nx-3) { //third-to-last column
                    temp[0] = data[row*Ny+col-3];
                    temp[1] = data[row*Ny+col-2];
                    temp[2] = data[row*Ny+col-1];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[row*Ny+col+1];
                    temp[5] = data[row*Ny+col+2];
                    temp[6] = data[row*Ny+0];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                }
                
            }
        }
        delete[] temp;
    };
    void derYFor(double* data, double* derivative) {
        double* temp = new double[7];
        for (int col = 0; col < Nx; col++) {
            for (int row = 0; row < Ny; row++) {
                if (3<=row && row<=Ny-4) {
                    temp[0] = data[(row-3)*Ny+col];
                    temp[1] = data[(row-2)*Ny+col];
                    temp[2] = data[(row-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //multiplied by zero
                    temp[4] = data[(row+1)*Ny+col];
                    temp[5] = data[(row+2)*Ny+col];
                    temp[6] = data[(row+3)*Ny+col];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (row == 0) { //first row
                    temp[0] = data[(Ny-3)*Ny+col];
                    temp[1] = data[(Ny-2)*Ny+col];
                    temp[2] = data[(Ny-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[(row+1)*Ny+col];
                    temp[5] = data[(row+2)*Ny+col];
                    temp[6] = data[(row+3)*Ny+col];

                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (row == 1) { //second row
                    temp[0] = data[(Ny-2)*Ny+col];
                    temp[1] = data[(Ny-1)*Ny+col];
                    temp[2] = data[(row-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[(row+1)*Ny+col];
                    temp[5] = data[(row+2)*Ny+col];
                    temp[6] = data[(row+3)*Ny+col];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (row == 2) { //third row
                    temp[0] = data[(Ny-1)*Ny+col];
                    temp[1] = data[(row-2)*Ny+col];
                    temp[2] = data[(row-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[(row+1)*Ny+col];
                    temp[5] = data[(row+2)*Ny+col];
                    temp[6] = data[(row+3)*Ny+col];

                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (row == Nx-1) { //last row
                    temp[0] = data[(row-3)*Ny+col];
                    temp[1] = data[(row-2)*Ny+col];
                    temp[2] = data[(row-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[0*Ny+col];
                    temp[5] = data[1*Ny+col];
                    temp[6] = data[2*Ny+col];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (row == Ny-2) { //second-to-last row
                    temp[0] = data[(row-3)*Ny+col];
                    temp[1] = data[(row-2)*Ny+col];
                    temp[2] = data[(row-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[(row+1)*Ny+col];
                    temp[5] = data[0*Ny+col];
                    temp[6] = data[1*Ny+col];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                } else if (row == Ny-3) { //third-to-last row
                    temp[0] = data[(row-3)*Ny+col];
                    temp[1] = data[(row-2)*Ny+col];
                    temp[2] = data[(row-1)*Ny+col];
                    temp[3] = data[row*Ny+col]; //0
                    temp[4] = data[(row+1)*Ny+col];
                    temp[5] = data[(row+2)*Ny+col];
                    temp[6] = data[0*Ny+col];
                    
                    derivative[row*Ny+col] = cblas_ddot(7,stencil,1,temp,1);
                }
                
            }
        }
        delete[] temp;
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
    int Nx = 100;
    int Ny = 100;
    int ic = 4;
    
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
    
    double dx = myInstance.dx;
    
    
    // output to file
    cout << endl << "Writing output to file: ";
    
    ofstream vOut("output.txt", ios::out | ios::trunc);
    vOut.precision(5);
    for (int yInd = 0;yInd<Ny;yInd++) {
        for (int xInd = 0;xInd<Nx;xInd++) {
        vOut << xInd*dx << "\t"
            << yInd*dx << "\t"
            << myInstance.yn[yInd*Ny+xInd] << "\t"
            << myInstance.yn[yInd*Ny+xInd + Nx*Ny] << "\t"
            << myInstance.yn[yInd*Ny+xInd + 2*Nx*Ny] << "\t" << endl;
        }
    }
    vOut.close();
    cout << "DONE" << endl;
}
