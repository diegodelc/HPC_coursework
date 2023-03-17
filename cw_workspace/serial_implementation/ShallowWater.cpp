#include <iostream>
using namespace std;
#include <cmath>
#include <cblas.h>

#include "ShallowWater.h"

ShallowWater::ShallowWater(double dt_in,double T_in,
                            int Nx_in,int Ny_in,
                            int ic_in) {
            cout << "Setting parameters: ";
            dt = dt_in;
            T = T_in;
            Nx = Nx_in;
            Ny = Ny_in;
            ic = ic_in;
            dx = 1;
            yn = new double[3*Nx*Ny];
            
            //u = new double[Nx*Ny];
            //v = new double[Nx*Ny];
            //h = new double[Nx*Ny];
            cout << "DONE" << endl;
        };
        
void ShallowWater::SetInitialConditions() {
            /*
            This function initialises the u, v and h arrays to the values defined in the assignement brief
            
                u and v are always zero, so they are initialised separately from h
                h is initialised as one of four intial conditions, specified via the parameter 'ic'
            */
            
            cout << "Setting initial conditions: ";
            
            
            for (int x = 0;x<Nx; x++) {
                for (int y = 0;y<Ny; y++) {
                    yn[x*Ny + y] = 0; //u
                    yn[Nx*Ny + x*Ny + y] = 0; //v
                }
            }
            

            //  This is slower than it has to be, it is evaluating the if Nx*Ny times, 
            //  would be better to evaluate if and then go into for
            //double xdouble = 0;
            //double ydouble;
            for (int x = 0;x<Nx; x++) {
                //ydouble = 0;
                for (int y = 0;y<Ny; y++) {
                    if (ic == 1) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond1((double)x,(double)y*dx);
                    } else if (ic == 2) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond2((double)x*dx,(double)y*dx);
                    } else if (ic == 3) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond3((double)x*dx,(double)y*dx);
                    } else if (ic == 4) {
                        yn[2*Nx*Ny + x*Ny + y] = initialCond4((double)x*dx,(double)y*dx);
                        
                    }
                    //ydouble++;
                }
                //xdouble++;
            }
            //cout << endl;
            cout << "DONE" << endl;
            
        };
        
void ShallowWater::TimeIntegrate() {
            cout << endl << "Starting time integration:" << endl;
            double* dudx = new double[Nx*Ny];
            double* dudy = new double[Nx*Ny];
            
            double* dvdx = new double[Nx*Ny];
            double* dvdy = new double[Nx*Ny];
            
            double* dhdx = new double[Nx*Ny];
            double* dhdy = new double[Nx*Ny];
            
            
            
            double* temp = new double[3*Nx*Ny];
            double* allKs = new double[3*Nx*Ny];
            // time propagation (for or while)
            //int counter = 0;
            for (double time = 0; time < T; time += dt) { //T and dt are double, so easier to make time double than try and cast them to int
                //This is for debugging, it prints everything
                /*
                if (counter == 50) {
                    cout << setw(15) << "h";
                    
                    cout << setw(15) << "temp";
                    
                    cout << setw(15) << "dudx";
                    cout << setw(15) << "dudy";
                    
                    cout << setw(15) << "dvdx";
                    cout << setw(15) << "dvdy";
                    
                    cout << setw(15) << "dhdx";
                    cout << setw(15) << "dhdy" << endl << endl;
                    
                    for (int i=0;i<Nx;i++) {
                        for (int j=0; j<Ny;j++) {
                            cout << setw(15) << yn[2*Nx*Ny + i*Ny + j]; //h
                            
                            cout << setw(15) << temp[Nx*Ny + i*Ny + j];
                            
                            cout << setw(15) << dudx[i*Ny + j];
                            cout << setw(15) << dudy[i*Ny + j];
                            
                            cout << setw(15) << dvdx[i*Ny + j];
                            cout << setw(15) << dvdy[i*Ny + j];
                            
                            cout << setw(15) << dhdx[i*Ny + j];
                            cout << setw(15) << dhdy[i*Ny + j] << endl;
                            
                        
                        }
                    }
                }
                */
                //k1 = calcF(yn)
                calcF(  yn,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);                
                for (int i = 0;i<3*Nx*Ny;i++) {
                    allKs[i] = temp[i];
                }
                

                //k2 = calcF(yn + dt*k1/2);
                for (int i = 0;i<3*Nx*Ny;i++) {
                    temp[i] = yn[i] + (dt/2) * temp[i];
                }
                calcF(  temp,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);
                for (int i = 0;i<3*Nx*Ny;i++) {
                    allKs[i] += 2*temp[i];
                }
                
                
                //k3 = calcF(yn + dt*k2/2);
                for (int i = 0;i<3*Nx*Ny;i++) {
                    temp[i] = yn[i] + (dt/2)*temp[i];
                }
                calcF(  temp,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);
                        
                
                for (int i = 0;i<3*Nx*Ny;i++) {
                    allKs[i] += 2*temp[i];
                }
                
                //k4 = calcF(yn + dt*k3);
                for (int i = 0;i<3*Nx*Ny;i++) {
                    temp[i] = yn[i] + dt*temp[i];
                }
                calcF(  temp,
                        dudx,dudy,
                        dvdx,dvdy,
                        dhdx,dhdy,
                        temp);
                
                for (int i = 0;i<3*Nx*Ny;i++) {
                    allKs[i] += temp[i];
                }
                
                //yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt;
                for (int i = 0;i<3*Nx*Ny;i++) {
                    yn[i] += (dt/6)*allKs[i];
                }                
                
            };
            
            cout << "\tFinished iteration " << T/dt << "/" << T/dt << endl;

            delete[] dudx;
            delete[] dudy;
            
            delete[] dvdx;
            delete[] dvdy;
            
            delete[] dhdx;
            delete[] dhdy;
            
            delete[] allKs;
            delete[] temp;
        };
        

//double stencil[7]; //declared in ShallowWater.h
double ShallowWater::initialCond1(double x,double y) {
        return 10 + exp(-(x-50)*(x-50)/25);
        //return 10 + exp(-(x-5)*(x-5)/2.5); //for a 10 by 10 grid
    };
double ShallowWater::initialCond2(double x,double y) {
        return 10 + exp(-(y-50)*(y-50)/25);
    };
double ShallowWater::initialCond3(double x,double y) {
        return 10 + exp(-((x-50)*(x-50) + (y-50)*(y-50))/25);
    };
double ShallowWater::initialCond4(double x,double y) {
        return 10 + exp(-((x-25)*(x-25) + (y-25)*(y-25))/25) + exp(-((x-75)*(x-75)+ (y-75)*(y-75) )/25);
    };
    
void ShallowWater::multVectByConst(const int& n,double* vect1,const double& constVal,double* ans,const char& HOW) {
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
    
void ShallowWater::addTwoVectors(const int& n,double* vect1,double* vect2,const double& sign,double* ans) {
        for (int i = 0; i<n; i++) {
            ans[i] = vect1[i] + sign * vect2[i];
        }
    };
    
void ShallowWater::multVectByVect(const int& n,double* vect1,double* vect2,const double& sign,double* ans,const char& HOW) {
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
    
void ShallowWater::calcF( double* yn,
                double* dudx,double* dudy,
                double* dvdx,double* dvdy,
                double* dhdx,double* dhdy,
                double* f) {
        
        //cblas_dcopy(Nx*Ny,yn,1,u,1);
        //cblas_dcopy(Nx*Ny,yn+Nx*Ny,1,v,1);
        //cblas_dcopy(Nx*Ny,yn+2*Nx*Ny,1,h,1);
        
        
        //NEED TO CHANGE NAMES TO U_POS SO THEY DON'T OVERWRITE THE ATTRIBUTES
        //ENSURE WE ARE NOT OVERWRITING YN, WE NEED THIS INTACT FOR EACH K EVALUATION
        //double* u_pos = yn;
        //double* v_pos = yn+Nx*Ny;
        //double* h_pos = yn+2*Nx*Ny;
        
        
        derXFor(yn,dudx);
        derYFor(yn,dudy);
        
        derXFor(yn+Nx*Ny,dvdx);
        derYFor(yn+Nx*Ny,dvdy);
        
        derXFor(yn+2*Nx*Ny,dhdx);
        derYFor(yn+2*Nx*Ny,dhdy);
        
        //f1 = - (g*dhdx + u.*dudx) - (v.*dudy);
        for (int i = 0; i<Nx*Ny; i++) {
            //f[i] = -9.81*dhdx[i] - u_pos[i]*dudx[i] - v_pos[i]*dudy[i];
            f[i] = -9.81*dhdx[i] - yn[i]*dudx[i] - yn[Nx*Ny + i]*dudy[i];
        }
        //multVectByConst(Nx*Ny,dhdx,-9.81,f,'R');
        //multVectByVect(Nx*Ny,dudx,u_pos,-1,f,'A');
        //multVectByVect(Nx*Ny,dudy,v_pos,-1,f,'A');
        
        //f2 = - (u.*dvdx) - (g*dhdy + v.*dvdy);
        for (int i = 0; i<Nx*Ny; i++) {
            //f[Nx*Ny + i] = - u_pos[i]*dvdx[i] - 9.81*dhdy[i] - v_pos[i]*dvdy[i];
            f[Nx*Ny + i] = - yn[i]*dvdx[i] - 9.81*dhdy[i] - yn[Nx*Ny + i]*dvdy[i];
        }
        //multVectByVect(Nx*Ny,u_pos,dvdx,-1,f+Nx*Ny,'R');
        //multVectByConst(Nx*Ny,dhdy,-9.81,f+Nx*Ny,'A');
        //multVectByVect(Nx*Ny,v_pos,dvdy,-1,f+Nx*Ny,'A');
        
        //f3 = - (u.*dhdx + h.*dudx) - (v.*dhdy + h.*dvdy);
        for (int i = 0; i<Nx*Ny; i++) {
            //f[2*Nx*Ny + i] = - u_pos[i]*dhdx[i] - h_pos[i]*dudx[i] - v_pos[i]*dhdy[i] - h[i]*dvdy[i];
            f[2*Nx*Ny + i] = - yn[i]*dhdx[i] - yn[2*Nx*Ny + i]*dudx[i] - yn[Nx*Ny + i]*dhdy[i] - yn[2*Nx*Ny + i]*dvdy[i];
        }
        //multVectByVect(Nx*Ny,u_pos,dhdx,-1,f+2*Nx*Ny,'R');
        //multVectByVect(Nx*Ny,h_pos,dudx,-1,f+2*Nx*Ny,'A');
        //multVectByVect(Nx*Ny,v_pos,dhdy,-1,f+2*Nx*Ny,'A');
        //multVectByVect(Nx*Ny,h_pos,dvdy,-1,f+2*Nx*Ny,'A');
    };
    
void ShallowWater::derXFor(double* data, double* derivative) {
        double* tempDer = new double[7];
        for (int xcol = 0; xcol < Nx; xcol++) {
            for (int yrow = 0; yrow < Ny; yrow++) {
                if (3<=xcol && xcol<=Nx-4) {
                    tempDer[0] = data[(xcol-3)*Ny+yrow];
                    tempDer[1] = data[(xcol-2)*Ny+yrow];
                    tempDer[2] = data[(xcol-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(xcol+1)*Ny+yrow];
                    tempDer[5] = data[(xcol+2)*Ny+yrow];
                    tempDer[6] = data[(xcol+3)*Ny+yrow];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (xcol == 0) { //first column                    
                    tempDer[0] = data[(Nx-3)*Ny+yrow];
                    tempDer[1] = data[(Nx-2)*Ny+yrow];
                    tempDer[2] = data[(Nx-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(xcol+1)*Ny+yrow];
                    tempDer[5] = data[(xcol+2)*Ny+yrow];
                    tempDer[6] = data[(xcol+3)*Ny+yrow];

                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (xcol == 1) { //second column                    
                    tempDer[0] = data[(Nx-2)*Ny+yrow];
                    tempDer[1] = data[(Nx-1)*Ny+yrow];
                    tempDer[2] = data[(xcol-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(xcol+1)*Ny+yrow];
                    tempDer[5] = data[(xcol+2)*Ny+yrow];
                    tempDer[6] = data[(xcol+3)*Ny+yrow];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (xcol == 2) { //third column                    
                    tempDer[0] = data[(Nx-1)*Ny+yrow];
                    tempDer[1] = data[(xcol-2)*Ny+yrow];
                    tempDer[2] = data[(xcol-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(xcol+1)*Ny+yrow];
                    tempDer[5] = data[(xcol+2)*Ny+yrow];
                    tempDer[6] = data[(xcol+3)*Ny+yrow];

                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (xcol == Nx-1) { //last column                    
                    tempDer[0] = data[(xcol-3)*Ny+yrow];
                    tempDer[1] = data[(xcol-2)*Ny+yrow];
                    tempDer[2] = data[(xcol-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(0)*Ny+yrow];
                    tempDer[5] = data[(1)*Ny+yrow];
                    tempDer[6] = data[(2)*Ny+yrow];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (xcol == Nx-2) { //second-to-last column                    
                    tempDer[0] = data[(xcol-3)*Ny+yrow];
                    tempDer[1] = data[(xcol-2)*Ny+yrow];
                    tempDer[2] = data[(xcol-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(xcol+1)*Ny+yrow];
                    tempDer[5] = data[(0)*Ny+yrow];
                    tempDer[6] = data[(1)*Ny+yrow];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (xcol == Nx-3) { //third-to-last column
                    tempDer[0] = data[(xcol-3)*Ny+yrow];
                    tempDer[1] = data[(xcol-2)*Ny+yrow];
                    tempDer[2] = data[(xcol-1)*Ny+yrow];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[(xcol+1)*Ny+yrow];
                    tempDer[5] = data[(xcol+2)*Ny+yrow];
                    tempDer[6] = data[(0)*Ny+yrow];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                }
                
            }
        }
        delete[] tempDer;
    };
void ShallowWater::derYFor(double* data, double* derivative) {
        double* tempDer = new double[7];
        for (int xcol = 0; xcol < Nx; xcol++) {
            for (int yrow = 0; yrow < Ny; yrow++) {
                if (3<=yrow && yrow<=Ny-4) {
                    tempDer[0] = data[xcol*Ny+yrow-3];
                    tempDer[1] = data[xcol*Ny+yrow-2];
                    tempDer[2] = data[xcol*Ny+yrow-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+yrow+1];
                    tempDer[5] = data[xcol*Ny+yrow+2];
                    tempDer[6] = data[xcol*Ny+yrow+3];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (yrow == 0) { //first row                    
                    tempDer[0] = data[xcol*Ny+Ny-3];
                    tempDer[1] = data[xcol*Ny+Ny-2];
                    tempDer[2] = data[xcol*Ny+Ny-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+yrow+1];
                    tempDer[5] = data[xcol*Ny+yrow+2];
                    tempDer[6] = data[xcol*Ny+yrow+3];

                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (yrow == 1) { //second row
                    tempDer[0] = data[xcol*Ny+Ny-2];
                    tempDer[1] = data[xcol*Ny+Ny-1];
                    tempDer[2] = data[xcol*Ny+yrow-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+yrow+1];
                    tempDer[5] = data[xcol*Ny+yrow+2];
                    tempDer[6] = data[xcol*Ny+yrow+3];

                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (yrow == 2) { //third row
                    tempDer[0] = data[xcol*Ny+Ny-1];
                    tempDer[1] = data[xcol*Ny+yrow-2];
                    tempDer[2] = data[xcol*Ny+yrow-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+yrow+1];
                    tempDer[5] = data[xcol*Ny+yrow+2];
                    tempDer[6] = data[xcol*Ny+yrow+3];


                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (yrow == Nx-1) { //last row
                    tempDer[0] = data[xcol*Ny+yrow-3];
                    tempDer[1] = data[xcol*Ny+yrow-2];
                    tempDer[2] = data[xcol*Ny+yrow-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+0];
                    tempDer[5] = data[xcol*Ny+1];
                    tempDer[6] = data[xcol*Ny+2];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (yrow == Ny-2) { //second-to-last row                    
                    tempDer[0] = data[xcol*Ny+yrow-3];
                    tempDer[1] = data[xcol*Ny+yrow-2];
                    tempDer[2] = data[xcol*Ny+yrow-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+yrow+1];
                    tempDer[5] = data[xcol*Ny+0];
                    tempDer[6] = data[xcol*Ny+1];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                } else if (yrow == Ny-3) { //third-to-last row
                    tempDer[0] = data[xcol*Ny+yrow-3];
                    tempDer[1] = data[xcol*Ny+yrow-2];
                    tempDer[2] = data[xcol*Ny+yrow-1];
                    tempDer[3] = data[xcol*Ny+yrow]; //multiplied by zero
                    tempDer[4] = data[xcol*Ny+yrow+1];
                    tempDer[5] = data[xcol*Ny+yrow+2];
                    tempDer[6] = data[xcol*Ny+0];
                    
                    derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
                }
                
            }
        }
        delete[] tempDer;
    };
    
