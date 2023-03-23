#include <iostream>
using namespace std;
#include <cmath>
#include <cblas.h>
#include <iomanip>

#include "ShallowWater.h"

#include <omp.h>

ShallowWater::ShallowWater(double dt_in,double T_in,
                            int Nx_in,int Ny_in,
                            int ic_in,
                            int integrationType_in) {
            cout << endl << "Setting parameters: ";
            dt = dt_in;
            T = T_in;
            Nx = Nx_in;
            Ny = Ny_in;
            ic = ic_in;
            dx = 1;
            yn = new double[3*Nx*Ny];
            integrationType = integrationType_in;
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
            
            cout << endl << "Setting initial conditions: ";
            
            #pragma omp parallel for
            for (int x = 0;x<Nx; x++) {
                for (int y = 0;y<Ny; y++) {
                    yn[x*Ny + y] = 0; //u
                    yn[Nx*Ny + x*Ny + y] = 0; //v
                }
            }
            

            
            
            if (ic == 1) {
                #pragma omp parallel for
                for (int x = 0;x<Nx; x++) {
                    double xd = (double)x;
                    //ydouble = 0;
                    for (int y = 0;y<Ny; y++) {
                        yn[2*Nx*Ny + x*Ny + y] = 10 + exp(-(xd-50)*(xd-50)/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
            } else if (ic == 2) {
                #pragma omp parallel for
                for (int y = 0;y<Ny; y++) {
                    double yd = (double)y;
                    //ydouble = 0;
                    for (int x = 0;x<Nx; x++) {
                        yn[2*Nx*Ny + x*Ny + y] = 10 + exp(-(yd-50)*(yd-50)/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
            } else if (ic == 3) {
                #pragma omp parallel for
                for (int x = 0;x<Nx; x++) {
                    double xd = (double)x;
                    //ydouble = 0;
                    for (int y = 0;y<Ny; y++) {
                        double yd = (double)y;
                        yn[2*Nx*Ny + x*Ny + y] = 10 + exp(-((xd-50)*(xd-50) + (yd-50)*(yd-50))/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
            } else if (ic == 4) {
                #pragma omp parallel for
                for (int x = 0;x<Nx; x++) {
                    double xd = (double)x;
                //ydouble = 0;
                    for (int y = 0;y<Ny; y++) {
                        double yd = (double)y;
                        yn[2*Nx*Ny + x*Ny + y] = 10 + exp(-((xd-25)*(xd-25) + (yd-25)*(yd-25))/25) + exp(-((xd-75)*(xd-75)+ (yd-75)*(yd-75) )/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
        }
            
            //cout << endl;
            cout << "DONE" << endl;
            
        };
        
void ShallowWater::TimeIntFor() {
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
        
        //k1 = calcF(yn)
        
        calcFFor(  yn,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                temp);                
        
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] = temp[i];
        //}
        

        //k2 = calcF(yn + dt*k1/2);
        //for (i = 0;i<3*Nx*Ny;i++) {
            temp[i] = yn[i] + (dt/2) * temp[i];
        }
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                temp);
                
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] += 2*temp[i];
        //}
        
        
        //k3 = calcF(yn + dt*k2/2);
        
        //for (i = 0;i<3*Nx*Ny;i++) {
            temp[i] = yn[i] + (dt/2)*temp[i];
        }
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                temp);
                
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] += 2*temp[i];
        //}
        
        //k4 = calcF(yn + dt*k3);
        
        //for (i = 0;i<3*Nx*Ny;i++) {
            temp[i] = yn[i] + dt*temp[i];
        }
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                temp);
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] += temp[i];
        //}
        
        //yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt;

        //for (i = 0;i<3*Nx*Ny;i++) {
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
void ShallowWater::calcFBLAS( double* yn,
                            double* dudx, double* dudy,
                            double* dvdx, double* dvdy,
                            double* dhdx, double* dhdy,
                            double* ddx,double* ddy,
                            double* A, double* B,
                            double* workspace,
                            double* F,
                            double* derXMat, double* derYMat,
                            double* vect,double* ans) {
    derXBlas(yn,dudx,derXMat,vect,ans);
    derYBlas(yn,dudy,derYMat,vect,ans);
    
    
    derXBlas(yn+Nx*Ny,dvdx,derXMat,vect,ans);
    derYBlas(yn+Nx*Ny,dvdy,derYMat,vect,ans);
    
    derXBlas(yn+2*Nx*Ny,dhdx,derXMat,vect,ans);
    derYBlas(yn+2*Nx*Ny,dhdy,derYMat,vect,ans);
    
    
    
    //stack derivatives in correct orders
    #pragma omp parallel for
    for (int i = 0; i<Nx*Ny; i++) {
        //stack x derivatives
        ddx[i*3] = dudx[i];
        ddx[i*3+1] = dvdx[i];
        ddx[i*3+2] = dhdx[i];
        
        //stack y derivatives
        ddy[i*3] = dudy[i];
        ddy[i*3+1] = dvdy[i];
        ddy[i*3+2] = dhdy[i];
    };
    
    //build A
    int kla = 2;
    int kua = 2;
    int lda = 1 + kla + kua;
    //double* A = new double[lda*Nx*Ny*3];
    #pragma omp parallel for
    for (int i = 0; i<Nx*Ny; i++) {
        //first column
        A[i*lda*3] = 0;
        A[i*lda*3+1] = 0;
        A[i*lda*3+2] = -yn[i];
        A[i*lda*3+3] = 0;
        A[i*lda*3+4] = -yn[2*Nx*Ny + i];
        
        //second column
        A[i*lda*3+5] = 0;
        A[i*lda*3+6] = 0;
        A[i*lda*3+7] = -yn[i];
        A[i*lda*3+8] = 0;
        A[i*lda*3+9] = 0;
        
        //third column
        A[i*lda*3+10] = -9.81;
        A[i*lda*3+11] = 0;
        A[i*lda*3+12] = -yn[i];
        A[i*lda*3+13] = 0;
        A[i*lda*3+14] = 0;
    }
    
    //build B
    int klb = 1;
    int kub = 1;
    int ldb = 1 + klb + kub;
    //double* B = new double[ldb*Nx*Ny*3];
    #pragma omp parallel for
    for (int i = 0; i<Nx*Ny; i++) {
        //first column
        B[i*ldb*3] = 0;
        B[i*ldb*3+1] = -yn[Nx*Ny + i];
        B[i*ldb*3+2] = 0;
        
        //second column
        B[i*ldb*3+3] = 0;
        B[i*ldb*3+4] = -yn[Nx*Ny + i];
        B[i*ldb*3+5] = -yn[2*Nx*Ny + i];
        
        //third column
        B[i*ldb*3+6] = -9.81;
        B[i*ldb*3+7] = -yn[Nx*Ny + i];
        B[i*ldb*3+8] = 0;
        
    }
    
    //solve system
    //cout << "Before first blas" << endl;
    cblas_dgbmv(CblasColMajor,CblasNoTrans,3*Nx*Ny,3*Nx*Ny,kla,kua,1,A,lda,ddx,1,0,F,1);
    
    //cout << "Before second blas" << endl;
    cblas_dgbmv(CblasColMajor,CblasNoTrans,3*Nx*Ny,3*Nx*Ny,klb,kub,1,B,ldb,ddy,1,1,F,1);
    
    
    //reshape F back to [u;v;h] shape for output
    cblas_dcopy(3*Nx*Ny,F,1,workspace,1);
    #pragma omp parallel for
    for (int i = 0; i<Nx*Ny; i++) {
        //u
        F[i] = workspace[i*3];
        
        //v
        F[Nx*Ny + i] = workspace[i*3 + 1];
        
        //h
        F[2*Nx*Ny + i] = workspace[i*3 + 2];
    }
    
}

void ShallowWater::derXBlas(double* data, double* derivative, double* derMat,double* vect, double* ans) {
    //NOTE: data comes in column major
    
    
    //these are being declared a second time, find a better way (maybe make them attributes?)
    int paddedLen = Nx+6;
    int kl = 3;
    int ku = 3;
    int lda = 1 + kl + ku;
    for (int yrow = 0; yrow<Ny; yrow++) {
        //first three (padding)
        vect[0] = data[(Ny-3)*Ny + yrow];//pad with third to last
        vect[1] = data[(Ny-2)*Ny + yrow];//pad with second to last
        vect[2] = data[(Ny-1)*Ny + yrow];//pad with last
        
        //last three (padding)
        vect[paddedLen-3] = data[(0)*Ny + yrow]; //pad with first
        vect[paddedLen-2] = data[(1)*Ny + yrow]; //pad with 2nd
        vect[paddedLen-1] = data[(2)*Ny + yrow]; //pad with 3rd
        
        //populate rest of vector
        for (int xcol = 0; xcol<Nx;xcol++) {
            vect[3+xcol] = data[(xcol)*Ny + yrow];
        }
        
        //blas operation 
        cblas_dgbmv(CblasColMajor,CblasNoTrans,paddedLen,paddedLen,kl,ku,1,derMat,lda,vect,1,0,ans,1);
        
        //map results to derivatives array (output)
        for (int xcol = 0; xcol<Nx; xcol++) {
            derivative[xcol*Ny+yrow] = ans[3+xcol];
        }
        
        
    }
    
}
void ShallowWater::derYBlas(double* data, double* derivative, double* derMat,double* vect, double* ans) {
    //NOTE: data comes in column major
    
    
    //these are being declared a second time, find a better way (maybe make them attributes?)
    int paddedLen = Ny+6;
    int kl = 3;
    int ku = 3;
    int lda = 1 + kl + ku;
    
    for (int xcol = 0; xcol<Nx; xcol++) {
        
        //first three (padding)
        vect[0] = data[(xcol)*Ny + Ny-1 -2];//pad with third to last
        vect[1] = data[(xcol)*Ny + Ny-1 -1];//pad with second to last
        vect[2] = data[(xcol)*Ny + Ny-1];//pad with last
        
        //last three (padding)
        vect[paddedLen-3] = data[(xcol)*Ny]; //pad with first
        vect[paddedLen-2] = data[(xcol)*Ny + 1]; //pad with 2nd
        vect[paddedLen-1] = data[(xcol)*Ny + 2]; //pad with 3rd
        
        //populate rest of vector
        for (int yrow = 0; yrow<Ny;yrow++) {
            vect[3+yrow] = data[(xcol)*Ny + yrow];
        }
        
        //blas operation 
        cblas_dgbmv(CblasColMajor,CblasNoTrans,paddedLen,paddedLen,kl,ku,1,derMat,lda,vect,1,0,ans,1);
        
        //map results to derivatives array (output)
        for (int yrow = 0; yrow<Ny; yrow++) {
            derivative[xcol*Ny+yrow] = ans[3+yrow];
        }
        
        
    }
    
}

void ShallowWater::TimeIntBlas() {
    //cout << "Starting BLAS time integration" << endl;
    
    
    //calculate derivatives (these are being done with for loops, could also be calculated using blas routines)
    double* dudx = new double[Nx*Ny];
    double* dudy = new double[Nx*Ny];
    
    double* dvdx = new double[Nx*Ny];
    double* dvdy = new double[Nx*Ny];
    
    double* dhdx = new double[Nx*Ny];
    double* dhdy = new double[Nx*Ny];
    
    double* ddx = new double[3*Nx*Ny];
    double* ddy = new double[3*Nx*Ny];
    
    
    double* temp = new double[3*Nx*Ny];
    double* workspace = new double[3*Nx*Ny];
    double* allKs = new double[3*Nx*Ny];
    
    double* A = new double[5*3*Nx*Ny];
    double* B = new double[3*3*Nx*Ny];
    
    
    
    
    //these need to be large enough for both the X and Y derivative
    int paddedLen;
    
    if (Ny>Nx) {
        paddedLen = Ny+6;
    } else {
        paddedLen = Nx+6;
    }
    double* vect = new double[paddedLen];
    double* ans = new double[paddedLen];
    //creating differentiation matrix for both X and Y
    int kl = 3;
    int ku = 3;
    int lda = 1 + ku + kl;
    
    int paddedLenX = Nx + 6;
    double* derXMat = new double[lda*paddedLenX]; //dimensions corresponding to largest grid dimension
    #pragma omp parallel for
    for (int i = 0; i<paddedLenX;i++) {
        derXMat[i*lda]   = -0.0167;     //c7
        derXMat[i*lda+1] =  0.1500;     //c6
        derXMat[i*lda+2] = -0.7500;     //c5
        derXMat[i*lda+3] =  0;          //c4
        derXMat[i*lda+4] =  0.7500;     //c3
        derXMat[i*lda+5] = -0.1500;     //c2
        derXMat[i*lda+6] =  0.0167;     //c1
    }
    
    int paddedLenY = Ny + 6;
    double* derYMat = new double[lda*paddedLenY]; //dimensions corresponding to largest grid dimension
    #pragma omp parallel for
    for (int i = 0; i<paddedLenY;i++) {
        derYMat[i*lda]   = -0.0167;     //c7
        derYMat[i*lda+1] =  0.1500;     //c6
        derYMat[i*lda+2] = -0.7500;     //c5
        derYMat[i*lda+3] =  0;          //c4
        derYMat[i*lda+4] =  0.7500;     //c3
        derYMat[i*lda+5] = -0.1500;     //c2
        derYMat[i*lda+6] =  0.0167;     //c1
    }
    
    // time propagation
    for (double time = 0; time < T; time += dt) {
        
        
        
        //k1 = calcF(yn)
        calcFBLAS(  yn,
                    dudx, dudy,
                    dvdx, dvdy,
                    dhdx, dhdy,
                    ddx, ddy,
                    A, B,
                    workspace,
                    temp,
                    derXMat,derYMat,
                    vect,ans);           
        
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] = temp[i];
        
            //k2 = calcF(yn + dt*k1/2);
            temp[i] = yn[i] + (dt/2) * temp[i];
        }
        calcFBLAS(  temp,
                    dudx, dudy,
                    dvdx, dvdy,
                    dhdx, dhdy,
                    ddx, ddy,
                    A, B,
                    workspace,
                    temp,
                    derXMat,derYMat,
                    vect,ans); 
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] += 2*temp[i];
        
        
        
            //k3 = calcF(yn + dt*k2/2);
            temp[i] = yn[i] + (dt/2)*temp[i];
        }
        calcFBLAS(  temp,
                    dudx, dudy,
                    dvdx, dvdy,
                    dhdx, dhdy,
                    ddx, ddy,
                    A, B,
                    workspace,
                    temp,
                    derXMat,derYMat,
                    vect,ans);  
                
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] += 2*temp[i];
        
        
            //k4 = calcF(yn + dt*k3);
            temp[i] = yn[i] + dt*temp[i];
        }
        calcFBLAS(  temp,
                    dudx, dudy,
                    dvdx, dvdy,
                    dhdx, dhdy,
                    ddx, ddy,
                    A, B,
                    workspace,
                    temp,
                    derXMat,derYMat,
                    vect,ans); 
        #pragma omp parallel for
        for (int i = 0;i<3*Nx*Ny;i++) {
            allKs[i] += temp[i];
        
        
            //yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt;
            yn[i] += (dt/6)*allKs[i];
        }    
        
    }
    cout << "\tFinished iteration " << T/dt << "/" << T/dt << endl;
    
    //delete allocations
    delete[] dudx;
    delete[] dudy;
    
    delete[] dvdx;
    delete[] dvdy;
    
    delete[] dhdx;
    delete[] dhdy;
    
    delete[] ddx;
    delete[] ddy;
    
    delete[] temp;
    delete[] workspace;
    delete[] allKs;
    
    delete[] A;
    delete[] B;
    
    delete[] vect;
    delete[] ans;
    
    delete[] derXMat;
    delete[] derYMat;
};

void ShallowWater::TimeIntegrate() {
            cout << endl << "Starting time integration:" << endl;
            
            if (integrationType == 1) {
                cout << "\tFor loop implementation of time integration chosen" << endl;
                TimeIntFor();
            } else if (integrationType == 2) {
                cout << "\tBLAS implementation of time integration chosen" << endl;
                TimeIntBlas();
            } else {
                cout << "\tERROR: integration type <" << integrationType << "> not implemented" << endl;
                
            }
            
        };
        


    
void ShallowWater::calcFFor( double* yn,
                double* dudx,double* dudy,
                double* dvdx,double* dvdy,
                double* dhdx,double* dhdy,
                double* f) {
        
        
        //#pragma omp parallel
        //{
        //#pragma omp sections
        //{
            //#pragma omp section    
            derXFor(yn,dudx);
            //#pragma omp section    
            derYFor(yn,dudy);
            //#pragma omp section  
            derXFor(yn+Nx*Ny,dvdx);
            //#pragma omp section    
            derYFor(yn+Nx*Ny,dvdy);
            //#pragma omp section  
            derXFor(yn+2*Nx*Ny,dhdx);
            //#pragma omp section    
            derYFor(yn+2*Nx*Ny,dhdy);
        //}
        //}
        //f1 = - (g*dhdx + u.*dudx) - (v.*dudy);
        
        #pragma omp parallel for
        for (int i = 0; i<Nx*Ny; i++) {
            f[i] = -9.81*dhdx[i] - yn[i]*dudx[i] - yn[Nx*Ny + i]*dudy[i];
        
        
        
        //f2 = - (u.*dvdx) - (g*dhdy + v.*dvdy);
            f[Nx*Ny + i] = - yn[i]*dvdx[i] - 9.81*dhdy[i] - yn[Nx*Ny + i]*dvdy[i];
        
        
        
        //f3 = - (u.*dhdx + h.*dudx) - (v.*dhdy + h.*dvdy);
            f[2*Nx*Ny + i] = - yn[i]*dhdx[i] - yn[2*Nx*Ny + i]*dudx[i] - yn[Nx*Ny + i]*dhdy[i] - yn[2*Nx*Ny + i]*dvdy[i];
        }
        
    };
    
void ShallowWater::derXFor(const double* data, double* derivative) {
        #pragma omp parallel for
        for (int yrow = 0; yrow <Ny; yrow++) {
            //derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
            //0th
            derivative[yrow] =      data[(Nx-3)*Ny+yrow] * (-0.0167) +
                                    data[(Nx-2)*Ny+yrow] * 0.1500 +
                                    data[(Nx-1)*Ny+yrow] * (-0.7500) +
                                    data[(1)*Ny+yrow]    * 0.7500 +
                                    data[(2)*Ny+yrow]    * (-0.1500) +
                                    data[(3)*Ny+yrow]    * 0.0167;
            //1st
            derivative[1*Ny+yrow] = data[(Nx-2)*Ny+yrow] * (-0.0167) +
                                    data[(Nx-1)*Ny+yrow] * 0.1500 +
                                    data[(0)*Ny+yrow]    * (-0.7500) +
                                    data[(2)*Ny+yrow]    * 0.7500 +
                                    data[(3)*Ny+yrow]    * (-0.1500) +
                                    data[(4)*Ny+yrow]    * 0.0167;
            //2nd
            derivative[2*Ny+yrow] = data[(Nx-1)*Ny+yrow] * (-0.0167) +
                                    data[(0)*Ny+yrow]    * 0.1500 +
                                    data[(1)*Ny+yrow]    * (-0.7500) +
                                    data[(3)*Ny+yrow]    * 0.7500 +
                                    data[(4)*Ny+yrow]    * (-0.1500) +
                                    data[(5)*Ny+yrow]    * 0.0167;
            
            //last
            derivative[(Nx-1)*Ny+yrow] =   data[(Nx-4)*Ny+yrow] * (-0.0167) +
                                           data[(Nx-3)*Ny+yrow] * 0.1500 +
                                           data[(Nx-2)*Ny+yrow] * (-0.7500) +
                                           data[(0)*Ny+yrow]    * 0.7500 +
                                           data[(1)*Ny+yrow]    * (-0.1500) +
                                           data[(2)*Ny+yrow]    * 0.0167;
            //second to last
            derivative[(Nx-2)*Ny+yrow] =   data[(Nx-5)*Ny+yrow] * (-0.0167) +
                                           data[(Nx-4)*Ny+yrow] * 0.1500 +
                                           data[(Nx-3)*Ny+yrow] * (-0.7500) +
                                           data[(Nx-1)*Ny+yrow] * 0.7500 +
                                           data[(0)*Ny+yrow]    * (-0.1500) +
                                           data[(1)*Ny+yrow]    * 0.0167;
            //third to last
            derivative[(Nx-3)*Ny+yrow] =   data[(Nx-6)*Ny+yrow] * (-0.0167) +
                                           data[(Nx-5)*Ny+yrow] * 0.1500 +
                                           data[(Nx-4)*Ny+yrow] * (-0.7500) +
                                           data[(Nx-2)*Ny+yrow] * 0.7500 +
                                           data[(Nx-1)*Ny+yrow] * (-0.1500) +
                                           data[(0)*Ny+yrow]    * 0.0167;
            
            for (int xcol = 3; xcol < Nx-3; xcol++) {
                derivative[xcol*Ny+yrow] =     data[(xcol-3)*Ny+yrow] * (-0.0167) +
                                               data[(xcol-2)*Ny+yrow] * 0.1500 +
                                               data[(xcol-1)*Ny+yrow] * (-0.7500) +
                                               data[(xcol+1)*Ny+yrow] * 0.7500 +
                                               data[(xcol+2)*Ny+yrow] * (-0.1500) +
                                               data[(xcol+3)*Ny+yrow] * 0.0167;
            }                                                         

        }
        
    };
void ShallowWater::derYFor(const double* data, double* derivative) {
    #pragma omp parallel for
    for (int xcol = 0; xcol<Nx; xcol++) {
            //0th
            derivative[xcol*Ny] =  data[(xcol)*Ny+Ny-3] * (-0.0167) +
                                   data[(xcol)*Ny+Ny-2] * 0.1500 +
                                   data[(xcol)*Ny+Ny-1] * (-0.7500) +
                                   data[(xcol)*Ny+1] * 0.7500 +
                                   data[(xcol)*Ny+2] * (-0.1500) +
                                   data[(xcol)*Ny+3]    * 0.0167;
            //1st
            derivative[xcol*Ny+1] =    data[(xcol)*Ny+Ny-2] * (-0.0167) +
                                       data[(xcol)*Ny+Ny-1] * 0.1500 +
                                       data[(xcol)*Ny+0] * (-0.7500) +
                                       data[(xcol)*Ny+2] * 0.7500 +
                                       data[(xcol)*Ny+3] * (-0.1500) +
                                       data[(xcol)*Ny+4]    * 0.0167;
            //2nd
            derivative[xcol*Ny+2] =    data[(xcol)*Ny+Ny-1] * (-0.0167) +
                                       data[(xcol)*Ny+0] * 0.1500 +
                                       data[(xcol)*Ny+1] * (-0.7500) +
                                       data[(xcol)*Ny+3] * 0.7500 +
                                       data[(xcol)*Ny+4] * (-0.1500) +
                                       data[(xcol)*Ny+5]    * 0.0167;
            
            //last
            derivative[xcol*Ny+Ny-1] = data[(xcol)*Ny+Ny-4] * (-0.0167) +
                                       data[(xcol)*Ny+Ny-3] * 0.1500 +
                                       data[(xcol)*Ny+Ny-2] * (-0.7500) +
                                       data[(xcol)*Ny+0] * 0.7500 +
                                       data[(xcol)*Ny+1] * (-0.1500) +
                                       data[(xcol)*Ny+2]    * 0.0167;
            //second to last
            derivative[xcol*Ny+Ny-2] = data[(xcol)*Ny+Ny-5] * (-0.0167) +
                                       data[(xcol)*Ny+Ny-4] * 0.1500 +
                                       data[(xcol)*Ny+Ny-3] * (-0.7500) +
                                       data[(xcol)*Ny+Ny-1] * 0.7500 +
                                       data[(xcol)*Ny+0] * (-0.1500) +
                                       data[(xcol)*Ny+1]    * 0.0167;
            //third to last
            derivative[xcol*Ny+Ny-3] = data[(xcol)*Ny+Ny-6] * (-0.0167) +
                                       data[(xcol)*Ny+Ny-5] * 0.1500 +
                                       data[(xcol)*Ny+Ny-4] * (-0.7500) +
                                       data[(xcol)*Ny+Ny-2] * 0.7500 +
                                       data[(xcol)*Ny+Ny-1] * (-0.1500) +
                                       data[(xcol)*Ny+0]    * 0.0167;
            
            for (int yrow = 3;yrow<Ny-3;yrow++) {                
                derivative[xcol*Ny+yrow] = data[(xcol)*Ny+yrow-3] * (-0.0167) +
                                           data[(xcol)*Ny+yrow-2] * 0.1500 +
                                           data[(xcol)*Ny+yrow-1] * (-0.7500) +
                                           data[(xcol)*Ny+yrow+1] * 0.7500 +
                                           data[(xcol)*Ny+yrow+2] * (-0.1500) +
                                           data[(xcol)*Ny+yrow+3]    * 0.0167;

            }
        
        }
        
    };
    
