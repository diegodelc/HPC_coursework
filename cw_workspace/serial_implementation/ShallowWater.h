class ShallowWater 
{
public:
        double dt;
        double T;
        int Nx;
        int Ny;
        int ic;
        
        double dtover2;
        double dtover6;
        
        double* yn;
        
        int integrationType;
        
        double* u;
        double* v;
        double* h;
        
        double dx;
        
        int Nxy;
        int Nxy2;
        int Nxy3;
        
        ShallowWater(double dt_in,double T_in,
                            int Nx_in,int Ny_in,
                            int ic_in,
                            int integrationType_in);
        
        void SetInitialConditions();
        
        void TimeIntegrate();
        
private:

    
    double stencil[7] = {-0.0167, 0.1500, -0.7500, 0, 0.7500, -0.1500, 0.0167}; //somehow declaring this here makes the code much faster (global scope, cache?)
    
    void calcFFor( double* yn,
                double* dudx,double* dudy,
                double* dvdx,double* dvdy,
                double* dhdx,double* dhdy,
                double* f);
    
    void derXFor(const double* data, double* derivative);
    
    void derYFor(const double* data, double* derivative);
    
    void TimeIntFor();
    
    void TimeIntBlas();

    void calcFBLAS(double* yn,
                    double* dudx, double* dudy,
                    double* dvdx, double* dvdy,
                    double* dhdx, double* dhdy,
                    double* ddx,double* ddy,
                    double* A, double* B,
                    double* workspace,
                    double* F,
                    double* derMat,double* derYMat,
                    double* vect,double* ans);
    
    void derXBlas(double* data, double* derivative, double* derMat,double* vect,double* ans);
    
    void derYBlas(double* data, double* derivative, double* derMat,double* vect,double* ans);

                    
};