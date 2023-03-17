
class ShallowWater 
{
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
        
        double dx;
        
        ShallowWater(double dt_in,double T_in,
                            int Nx_in,int Ny_in,
                            int ic_in);
        
        void SetInitialConditions();
        
        void TimeIntegrate();
        
private:
    double stencil[7] = {-0.0167, 0.1500, -0.7500, 0, 0.7500, -0.1500, 0.0167}; //somehow declaring this here makes the code much faster
    double initialCond1(double x,double y);
    double initialCond2(double x,double y);
    double initialCond3(double x,double y);
    double initialCond4(double x,double y);
    
    void multVectByConst(const int& n,double* vect1,const double& constVal,double* ans,const char& HOW);
    
    void addTwoVectors(const int& n,double* vect1,double* vect2,const double& sign,double* ans);
    
    void multVectByVect(const int& n,double* vect1,double* vect2,const double& sign,double* ans,const char& HOW);
    
    void calcFFor( double* yn,
                double* dudx,double* dudy,
                double* dvdx,double* dvdy,
                double* dhdx,double* dhdy,
                double* f);
    
    void derXFor(double* data, double* derivative);
    
    void derYFor(double* data, double* derivative);
};