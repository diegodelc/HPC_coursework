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
                            int integrationType_in,
                            int verb_in) {
            verb = verb_in;
            if (verb == 1) {
                cout << endl << "Setting parameters: ";
            }
            dt = dt_in;
            dtover2 = dt/2;
            dtover6 = dt/6;
            T = T_in;
            Nx = Nx_in;
            Ny = Ny_in;
            Nxy = Nx*Ny;
            Nxy2 = 2*Nx*Ny;
            Nxy3 = 3*Nx*Ny;
            ic = ic_in;
            dx = 1;
            yn = new double[Nxy3];
            integrationType = integrationType_in;
            //u = new double[Nx*Ny];
            //v = new double[Nx*Ny];
            //h = new double[Nx*Ny];
            if (verb == 1) {
                cout << "DONE" << endl;
            }
            
        };
        
void ShallowWater::SetInitialConditions() {
            /*
            This function initialises the u, v and h arrays to the values defined in the assignement brief
            
                u and v are always zero, so they are initialised separately from h
                h is initialised as one of four intial conditions, specified via the parameter 'ic'
            */
            
            if (verb == 1) {
                cout << endl << "Setting initial conditions: ";
            }
            #pragma omp parallel for
            for (int x = 0;x<Nx; x++) {
                
                
                for (int y = 0;y<Ny; y++) {
                    yn[x*Ny + y] = 0; //u
                    yn[Nxy + x*Ny + y] = 0; //v
                }
            }
            

            double xd;
            double yd;
            
            if (ic == 1) {
                #pragma omp parallel for
                for (int x = 0;x<Nx; x++) {
                    xd = (double)x;
                    
                    //ydouble = 0;
                    for (int y = 0;y<Ny; y++) {
                        yn[Nxy2 + x*Ny + y] = 10 + exp(-(xd-50)*(xd-50)/25); //h
                        //ydouble++;
                    }
                    //xdouble++;
                }
            } else if (ic == 2) {
                #pragma omp parallel for
                for (int y = 0;y<Ny; y++) {
                    yd = (double)y;
                    //ydouble = 0;
                    for (int x = 0;x<Nx; x++) {
                        yn[Nxy2 + x*Ny + y] = 10 + exp(-(yd-50)*(yd-50)/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
            } else if (ic == 3) {
                #pragma omp parallel for
                for (int x = 0;x<Nx; x++) {
                    xd = (double)x;
                    
                    for (int y = 0;y<Ny; y++) {
                        yd = (double)y;
                        yn[Nxy2 + x*Ny + y] = 10 + exp(-(
                                                        (xd-50)*(xd-50) 
                                                      + (yd-50)*(yd-50)
                                                        )/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
            } else if (ic == 4) {
                #pragma omp parallel for
                for (int x = 0;x<Nx; x++) {
                    xd = (double)x;
                    
                //ydouble = 0;
                    for (int y = 0;y<Ny; y++) {
                        yd = (double)y;
                        yn[Nxy2 + x*Ny + y] = 10    + exp(-((xd-25)*(xd-25) + (yd-25)*(yd-25))/25) 
                                                    + exp(-((xd-75)*(xd-75) + (yd-75)*(yd-75) )/25);
                        //ydouble++;
                    }
                    //xdouble++;
                }
        }
            
            //cout << endl;
            if (verb == 1) {
                cout << "DONE" << endl;
            }
        };
        
void ShallowWater::TimeIntFor() {
    double* dudx = new double[Nxy];
    double* dudy = new double[Nxy];
    
    double* dvdx = new double[Nxy];
    double* dvdy = new double[Nxy];
    
    double* dhdx = new double[Nxy];
    double* dhdy = new double[Nxy];
    
    
    
    double* temp = new double[3*Nx*Ny];
    
    double* k1 = new double[Nxy3];
    double* k2 = new double[Nxy3];
    double* k3 = new double[Nxy3];
    double* k4 = new double[Nxy3];
    
    //double* allKs = new double[3*Nx*Ny];
    // time propagation (for or while)
    int counter = 0;
    
    for (double time = 0; time < T+dtover2; time += dt) { //T and dt are double, so easier to make time double than try and cast them to int
        //This is for debugging, it prints everything
        
        //k1 = calcF(yn)
        for (int i = 0;i<Nxy3;i++) {
            temp[i] = yn[i];
        }
        
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                k1);                
        
        #pragma omp parallel for
        for (int i = 0;i<Nxy3;i++) {
            //allKs[i] = temp[i];
        //}
        

        //k2 = calcF(yn + dt*k1/2);
        //for (i = 0;i<3*Nx*Ny;i++) {
            //temp[i] = yn[i] + (dt/2) * temp[i];
            temp[i] = yn[i] + dtover2 * k1[i];
        }
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                k2);
                
        #pragma omp parallel for
        for (int i = 0;i<Nxy3;i++) {
            //allKs[i] += 2*temp[i];
        //}
        
        
        //k3 = calcF(yn + dt*k2/2);
        
        //for (i = 0;i<3*Nx*Ny;i++) {
            temp[i] = yn[i] + dtover2*k2[i];
        }
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                k3);
                
        #pragma omp parallel for
        for (int i = 0;i<Nxy3;i++) {
            //allKs[i] += 2*temp[i];
        //}
        
        //k4 = calcF(yn + dt*k3);
        
        //for (i = 0;i<3*Nx*Ny;i++) {
            temp[i] = yn[i] + dt*k3[i];
        }
        calcFFor(  temp,
                dudx,dudy,
                dvdx,dvdy,
                dhdx,dhdy,
                k4);
        #pragma omp parallel for
        for (int i = 0;i<Nxy3;i++) {
            //allKs[i] += temp[i];
        //}
        
        //yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt;

        //for (i = 0;i<3*Nx*Ny;i++) {
            yn[i] += (dtover6)*(k1[i] +2*k2[i] + 2*k3[i] + k4[i]);
        }        
        counter ++;
        
    };
    
    if (verb == 1) {
        cout << "\tFinished iteration " << counter-1 << "/" << T/dt << endl; //subtract one since added counter at end of loop
        //cout << "\tcounter: " << counter << endl;
    }
    
    delete[] dudx;
    delete[] dudy;
    
    delete[] dvdx;
    delete[] dvdy;
    
    delete[] dhdx;
    delete[] dhdy;
    
    //delete[] allKs;
    //delete[] temp;
    
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
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
    
    
    derXBlas(yn+Nxy,dvdx,derXMat,vect,ans);
    derYBlas(yn+Nxy,dvdy,derYMat,vect,ans);
    
    derXBlas(yn+Nxy2,dhdx,derXMat,vect,ans);
    derYBlas(yn+Nxy2,dhdy,derYMat,vect,ans);
    
    
    //stack derivatives in correct orders for blas
    #pragma omp parallel for
    for (int i = 0; i<Nxy; i++) {
        int index = i*3;
        //stack x derivatives
        ddx[index] = dudx[i];
        ddx[index+1] = dvdx[i];
        ddx[index+2] = dhdx[i];
        
        //stack y derivatives
        ddy[index] = dudy[i];
        ddy[index+1] = dvdy[i];
        ddy[index+2] = dhdy[i];
    };
    
    //build A
    int kla = 2;
    int kua = 2;
    int lda = 1 + kla + kua;
    //double* A = new double[lda*Nx*Ny*3];
    #pragma omp parallel for
    for (int i = 0; i<Nxy; i++) {
        int index = i*lda*3;
        //first column
        A[index] = 0;
        A[index+1] = 0;
        A[index+2] = -yn[i];
        A[index+3] = 0;
        A[index+4] = -yn[2*Nx*Ny + i];
        
        //second column
        A[index+5] = 0;
        A[index+6] = 0;
        A[index+7] = -yn[i];
        A[index+8] = 0;
        A[index+9] = 0;
        
        //third column
        A[index+10] = -9.81;
        A[index+11] = 0;
        A[index+12] = -yn[i];
        A[index+13] = 0;
        A[index+14] = 0;
    }
    
    //build B
    int klb = 1;
    int kub = 1;
    int ldb = 1 + klb + kub;
    //double* B = new double[ldb*Nx*Ny*3];
    #pragma omp parallel for
    for (int i = 0; i<Nx*Ny; i++) {
        int index = i*ldb*3;
        //first column
        B[index] = 0;
        B[index+1] = -yn[Nx*Ny + i];
        B[index+2] = 0;
        
        //second column
        B[index+3] = 0;
        B[index+4] = -yn[Nx*Ny + i];
        B[index+5] = -yn[2*Nx*Ny + i];
        
        //third column
        B[index+6] = -9.81;
        B[index+7] = -yn[Nx*Ny + i];
        B[index+8] = 0;
        
    }
    
    //solve system
    //cout << "Before first blas" << endl;
    cblas_dgbmv(CblasColMajor,CblasNoTrans,3*Nx*Ny,3*Nx*Ny,kla,kua,1,A,lda,ddx,1,0,F,1);
    
    //cout << "Before second blas" << endl;
    cblas_dgbmv(CblasColMajor,CblasNoTrans,3*Nx*Ny,3*Nx*Ny,klb,kub,1,B,ldb,ddy,1,1,F,1);
    
    
    //reshape F back to [u;v;h] shape for output
    cblas_dcopy(Nxy3,F,1,workspace,1);
    #pragma omp parallel for
    for (int i = 0; i<Nxy; i++) {
        //u
        F[i] = workspace[i*3];
        
        //v
        F[Nxy + i] = workspace[i*3 + 1];
        
        //h
        F[Nxy2 + i] = workspace[i*3 + 2];
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
        vect[0] = data[(Nx-3)*Ny + yrow];//pad with third to last
        vect[1] = data[(Nx-2)*Ny + yrow];//pad with second to last
        vect[2] = data[(Nx-1)*Ny + yrow];//pad with last
        
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
    double* dudx = new double[Nxy];
    double* dudy = new double[Nxy];
    
    double* dvdx = new double[Nxy];
    double* dvdy = new double[Nxy];
    
    double* dhdx = new double[Nxy];
    double* dhdy = new double[Nxy];
    
    double* ddx = new double[Nxy3];
    double* ddy = new double[Nxy3];
    
    
    double* temp = new double[Nxy3];
    double* workspace = new double[Nxy3];
    double* allKs = new double[Nxy3];
    
    double* A = new double[5*Nxy3];
    double* B = new double[3*Nxy3];
    
    
    
    
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
        int index = i*lda;
        
        derXMat[index]   = -0.0167;     //7
        derXMat[index+1] =  0.1500;     //6
        derXMat[index+2] = -0.7500;     //5
        derXMat[index+3] =  0;          //4
        derXMat[index+4] =  0.7500;     //3
        derXMat[index+5] = -0.1500;     //2
        derXMat[index+6] =  0.0167;     //1
    }
    
    int paddedLenY = Ny + 6;
    double* derYMat = new double[lda*paddedLenY]; //dimensions corresponding to largest grid dimension
    #pragma omp parallel for
    for (int i = 0; i<paddedLenY;i++) {
        int index = i*lda;
        derYMat[index]   = -0.0167;     //7
        derYMat[index+1] =  0.1500;     //6
        derYMat[index+2] = -0.7500;     //5
        derYMat[index+3] =  0;          //4
        derYMat[index+4] =  0.7500;     //3
        derYMat[index+5] = -0.1500;     //2
        derYMat[index+6] =  0.0167;     //1
    }
    
    // time propagation
    int counter = 0;
    for (double time = 0; time < T + dtover2; time += dt) {
        
        
        
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
        for (int i = 0;i<Nxy3;i++) {
            allKs[i] = temp[i];
        
            //k2 = calcF(yn + dt*k1/2);
            temp[i] = yn[i] + (dtover2) * temp[i];
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
        for (int i = 0;i<Nxy3;i++) {
            allKs[i] += 2*temp[i];
        
        
        
            //k3 = calcF(yn + dt*k2/2);
            temp[i] = yn[i] + (dtover2)*temp[i];
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
        for (int i = 0;i<Nxy3;i++) {
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
        for (int i = 0;i<Nxy3;i++) {
            allKs[i] += temp[i];
        
        
            //yn = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)*dt;
            yn[i] += (dtover6)*allKs[i];
        }    
        counter++;
    }
    if (verb == 1) {
        cout << "\tFinished iteration " << counter-1 << "/" << T/dt << endl; //subtract one since added counter at end of loop
        //cout << "\tcounter: " << counter << endl;
    }
    
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
            if (verb == 1) {
                cout << endl << "Starting time integration:" << endl;
            }
            if (integrationType == 1) {
                if (verb == 1) {
                    cout << "\tFor loop implementation of time integration chosen" << endl;
                }
                TimeIntFor();
            } else if (integrationType == 2) {
                if (verb == 1) {
                    cout << "\tBLAS implementation of time integration chosen" << endl;
                }
                TimeIntBlas();
            } else {
                if (verb == 1) {
                    cout << "\tERROR: integration type <" << integrationType << "> not implemented" << endl;
                }
                
            }
            
        };
        


    
void ShallowWater::calcFFor( double* yn,
                double* dudx,double* dudy,
                double* dvdx,double* dvdy,
                double* dhdx,double* dhdy,
                double* f) {
        
                    
        derXFor(yn,dudx);
        derYFor(yn,dudy);
        
        derXFor(yn+Nxy,dvdx);        
        derYFor(yn+Nxy,dvdy);
        
        derXFor(yn+Nxy2,dhdx);
        derYFor(yn+Nxy2,dhdy);
        
        
        #pragma omp parallel for
        for (int i = 0; i<Nxy; i++) {

            //f1 = - (g*dhdx + u.*dudx) - (v.*dudy);
            f[i] = -9.81*dhdx[i] - yn[i]*dudx[i] - yn[Nxy + i]*dudy[i];
        
        
        
            //f2 = - (u.*dvdx) - (g*dhdy + v.*dvdy);
            f[Nxy + i] = - yn[i]*dvdx[i] - 9.81*dhdy[i] - yn[Nxy + i]*dvdy[i];
        
        
        
            //f3 = - (u.*dhdx + h.*dudx) - (v.*dhdy + h.*dvdy);
            f[Nxy2 + i] = - yn[i]*dhdx[i] - yn[Nxy2 + i]*dudx[i] - yn[Nxy + i]*dhdy[i] - yn[Nxy2 + i]*dvdy[i];
        }
        
    };
    
void ShallowWater::derXFor(const double* data, double* derivative) {
        #pragma omp parallel for
        for (int yrow = 0; yrow <Ny; yrow++) {
            //int temp;
            //derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,tempDer,1);
            //0th
            derivative[yrow] =      data[(Nx-3)*Ny+yrow] * (-0.0167) +               
                                    data[(Nx-2)*Ny+yrow] * ( 0.1500) +
                                    data[(Nx-1)*Ny+yrow] * (-0.7500) +
                                    //data[(0)*Ny+yrow] * 0( +     )
                                    data[(1)*Ny+yrow]    * ( 0.7500) +
                                    data[(2)*Ny+yrow]    * (-0.1500) +
                                    data[(3)*Ny+yrow]    * ( 0.0167);
            //1st                                                 
            derivative[1*Ny+yrow] = data[(Nx-2)*Ny+yrow] * (-0.0167) +
                                    data[(Nx-1)*Ny+yrow] * ( 0.1500) +
                                    data[(0)*Ny+yrow]    * (-0.7500) +
                                    //data[(1)*Ny+yrow] * 0( +     )
                                    data[(2)*Ny+yrow]    * ( 0.7500) +
                                    data[(3)*Ny+yrow]    * (-0.1500) +
                                    data[(4)*Ny+yrow]    * ( 0.0167);
            //2nd
            derivative[2*Ny+yrow] = data[(Nx-1)*Ny+yrow] * (-0.0167) +
                                    data[(0)*Ny+yrow]    * ( 0.1500) +
                                    data[(1)*Ny+yrow]    * (-0.7500) +
                                    //data[(2)*Ny+yrow] * 0( +     )
                                    data[(3)*Ny+yrow]    * ( 0.7500) +
                                    data[(4)*Ny+yrow]    * (-0.1500) +
                                    data[(5)*Ny+yrow]    * ( 0.0167);
            
            //last
            derivative[(Nx-1)*Ny+yrow] =   data[(Nx-4)*Ny+yrow] * (-0.0167) +
                                           data[(Nx-3)*Ny+yrow] * ( 0.1500) +
                                           data[(Nx-2)*Ny+yrow] * (-0.7500) +
                                           //data[(Nx-1)*Ny+yrow] (* 0 +  )
                                           data[(0)*Ny+yrow]    * ( 0.7500) +
                                           data[(1)*Ny+yrow]    * (-0.1500) +
                                           data[(2)*Ny+yrow]    * ( 0.0167);
            //second to last
            derivative[(Nx-2)*Ny+yrow] =   data[(Nx-5)*Ny+yrow] * (-0.0167) +
                                           data[(Nx-4)*Ny+yrow] * ( 0.1500) +
                                           data[(Nx-3)*Ny+yrow] * (-0.7500) +
                                           //data[(Nx-2)*Ny+yrow] (* 0 +  )
                                           data[(Nx-1)*Ny+yrow] * ( 0.7500) +
                                           data[(0)*Ny+yrow]    * (-0.1500) +
                                           data[(1)*Ny+yrow]    * ( 0.0167);
            //third to last
            derivative[(Nx-3)*Ny+yrow] =   data[(Nx-6)*Ny+yrow] * (-0.0167) +
                                           data[(Nx-5)*Ny+yrow] * ( 0.1500) +
                                           data[(Nx-4)*Ny+yrow] * (-0.7500) +
                                           //data[(Nx-3)*Ny+yrow] (* 0 +  )
                                           data[(Nx-2)*Ny+yrow] * ( 0.7500) +
                                           data[(Nx-1)*Ny+yrow] * (-0.1500) +
                                           data[(0)*Ny+yrow]    * ( 0.0167);
            
            for (int xcol = 3; xcol < Nx-3; xcol++) {
                /*
                temp = (xcol)*Ny+yrow;
                derivative[xcol*Ny+yrow] =     data[temp-Ny-Ny-Ny] * (-0.0167) +
                                               data[temp-Ny-Ny] * 0.1500 +
                                               data[temp-Ny] * (-0.7500) +
                                               //data[temp] * (0) +
                                               data[temp+Ny] * 0.7500 +
                                               data[temp+Ny+Ny] * (-0.1500) +
                                               data[temp+Ny+Ny+Ny] * 0.0167;
                */
                derivative[xcol*Ny+yrow] =     data[(xcol-3)*Ny+yrow] * (-0.0167) +
                                               data[(xcol-2)*Ny+yrow] * ( 0.1500) +
                                               data[(xcol-1)*Ny+yrow] * (-0.7500) +
                                               //data[(xcol)*Ny+yrow] * (( 0) + )
                                               data[(xcol+1)*Ny+yrow] * ( 0.7500) +
                                               data[(xcol+2)*Ny+yrow] * (-0.1500) +
                                               data[(xcol+3)*Ny+yrow] * ( 0.0167);
                
            }                                                         

        }
        
    };
void ShallowWater::derYFor(const double* data, double* derivative) {
    #pragma omp parallel for
    for (int xcol = 0; xcol<Nx; xcol++) {
            int index = (xcol)*Ny;
            //0th
            derivative[index] =    data[index+Ny-3] * (-0.0167) +
                                   data[index+Ny-2] * ( 0.1500) +
                                   data[index+Ny-1] * (-0.7500) +
                                   //data[(xcol)*Ny+0]( * 0 + )
                                   data[index+1]    * ( 0.7500) +
                                   data[index+2]    * (-0.1500) +
                                   data[index+3]    * ( 0.0167);
            //1st
            derivative[index+1] =      data[index+Ny-2] * (-0.0167) +
                                       data[index+Ny-1] * ( 0.1500) +
                                       data[index+0]    * (-0.7500) +
                                       //data[(xcol)*Ny+1]( * 0 + )
                                       data[index+2]    * ( 0.7500) +
                                       data[index+3]    * (-0.1500) +
                                       data[index+4]    * ( 0.0167);
            //2nd
            derivative[index+2] =      data[index+Ny-1]  *( -0.0167) +
                                       data[index+0]     *(  0.1500) +
                                       data[index+1]     *( -0.7500) +
                                       //data[(xcol)*Ny+2]( * 0 +  )
                                       data[index+3]     *(  0.7500) +
                                       data[index+4]     *( -0.1500) +
                                       data[index+5]     *(  0.0167);
            //last                                        (        )
            derivative[index+Ny-1] =   data[index+Ny-4] * (-0.0167 )+
                                       data[index+Ny-3] * ( 0.1500 )+
                                       data[index+Ny-2] * (-0.7500 )+
                                       //data[(xcol)*Ny+Ny(-1] * 0 )+
                                       data[index+0]    * ( 0.7500 )+
                                       data[index+1]    * (-0.1500 )+
                                       data[index+2]    * ( 0.0167 );
            //second to last                              (        )
            derivative[index+Ny-2] =   data[index+Ny-5] * (-0.0167 )+
                                       data[index+Ny-4] * ( 0.1500 )+
                                       data[index+Ny-3] * (-0.7500 )+
                                       //data[(xcol)*Ny+Ny(-2] * 0 )+
                                       data[index+Ny-1] * ( 0.7500 )+
                                       data[index+0]    * (-0.1500 )+
                                       data[index+1]    * ( 0.0167 );
            //third to last                               (        )
            derivative[index+Ny-3] =   data[index+Ny-6]  *( -0.0167) +
                                       data[index+Ny-5]  *(  0.1500) +
                                       data[index+Ny-4]  *( -0.7500) +
                                       //data[(xcol)*Ny+  (        )
                                       data[index+Ny-2]  *(  0.7500) +
                                       data[index+Ny-1]  *( -0.1500) +
                                       data[index+0]     *(  0.0167);
            
            for (int yrow = 3;yrow<Ny-3;yrow++) {                
                
                derivative[index+yrow] =   data[index+yrow-3] * (-0.0167) +
                                           data[index+yrow-2] * ( 0.1500) +
                                           data[index+yrow-1] * (-0.7500) +
                                           data[index+yrow+1] * ( 0.7500) + 
                                           data[index+yrow+2] * (-0.1500) +
                                           data[index+yrow+3] * ( 0.0167);  
                                                                
                //derivative[xcol*Ny+yrow] = cblas_ddot(7,stencil,1,data+(xcol)*Ny+yrow-3,1);

            }
        
        }
        
    }
    
