//implementation file(.cpp file)contains implementation of function and class methods
#include "domain.h"
#include <vector>
#include <cstring>

//instances of class Rnd.
Rnd rnd;

//domain constructor
Domain::Domain(double x0, double dx, int ni)
    : x0(x0), dx(dx), ni(ni)
{
    // Allocate memory for phi, rho, and ef
    phi = new double[ni];
    rho = new double[ni];
    ef = new double[ni];
    xL = (ni-1)*dx;

    // Initialize arrays to zero using memset
    memset(phi, 0, ni * sizeof(double));
    memset(rho, 0, ni * sizeof(double));
    memset(ef, 0, ni * sizeof(double));
}  

// dealloacate memory when out of  scope.
Domain::~Domain()
{
    // Deallocate memory
    //std::cout<<"domain destructor called"<<std::endl;
    delete[] phi;
    delete[] rho;
    delete[] ef;
}


void Domain::display()
{
    display::print("cell spacing :",dx);
    display::print("system lenght :",xL);
    display::print("grid points :",ni);
    display::print("debye lenght :",LD);
    display::print("plasma frequency :",wp);
    display::print("Time step:",DT);
    display::print("boundary condition :",bc);
    display::print("plasma density :",density);
    display::print("electron thermal velocity  :",LD*wp,"m/s");
    display::print("beam streaming velocity  :",v_b*LD*wp,"m/s");
    display::print("gamma factor :",sqrt(1-(v_b*LD*wp/(3e8))*(v_b*LD*wp/(3e8))));
}

//set normalized parameter.
void Domain:: set_normparam(double LD, double wp)
{
    this->LD = LD;
    this->wp = wp;
    
}

//set time.
void Domain::set_time(double DT, int NUM_TS, int write_interval)
{
    this->DT = DT;
    this->NUM_TS = NUM_TS;
    this->write_interval = write_interval;
}

// set differetnt simulation parameter.
void Domain::set_simparam(double tempE, double tempI, double tempN, double tempB, double density, double v_e, double v_i, double v_n, double v_b, double alpha, int species_no)
{
    this->tempE = tempE;
    this->tempI = tempI;
    this->tempN = tempN;
    this->tempB = tempB;
    this->density = density;
    this->v_i = v_i;
    this->v_e = v_e;
    this->v_n = v_n;
    this->v_b = v_b;
    this->alpha = alpha;
    this->species_no = species_no;
}


//compute normalized charge density (n_e/N) at grid points.
void Domain::ComputeRho(vector<Species> &species)// "species" is conatiner(vector here) holding instances(many instances) of class "Species"
{
    //":" means range based for loop. here we are iterating over the elements contained in
    // "species" container using a variable sp(adress) and Species is kind of data type(we can also use auto if we we donot know data type )
    // here data type is complex not like int,float char. it is "Species".
    //we use address because we want to amke changes in the original list and donot want to  make a copy.
    //for(Species &sp : species)
    
    // Iterate over each grid point
    for (int i = 0; i < ni; i++) 
    {
        // Initialize density at current grid point to zero
        rho[i] = 0.0;
        // Iterate over each species in the vector
        for (Species &sp : species) 
        {
            // Add contribution of current species to rho at current grid point
            rho[i] += (1 / Const::QE) * (sp.charge * sp.den[i]);
        }
    }    
}

void Domain::ComputeRho1(Species &species1,Species &species2)
{
    for(int i = 0; i < ni; i++)
    {
        rho[i] = (1/Const::QE)*(species1.charge*species1.den[i] + species2.charge*species2.den[i]);
    }
}

//scatter densities at grid points
void Domain::Scatter(double lc, double value, double *field)
{
    int i = (int)lc;
    //cout<<ni<<endl;
    double di = lc-i;
    field[i] += value*(1-di);
    field[i+1] += value*(di);
}


//Interpolate density and eletric field value at grid points
double Domain::Interpolate(double lc, double *field)
{
    int i=(int)lc;
    double di = lc-i;
    double val = field[i]*(1-di) + field[i+1]*(di);
    return val;
}

//this function doesmot needed when the code is normalized
double Domain:: XtoL(double pos)
{
    double li = (pos-x0)/dx;
    //cout<<li<<endl;
    return li;
}


//function to print 1d array/matrix
void Domain::printmat(double *A) 
{
    for (int i = 0; i < ni; i++) 
    {
        std::cout << A[i] << " ";
    }
    std::cout << std::endl;
}

double Domain::ComputePE(Species &species)
{
    double pe = 0;
    double eft = 0;
    for (int i= 0; i < ni; i++)
    {
        eft += ef[i];
    }
    // un-normalize this electric field
    eft = eft*(Const::ME*wp*wp*LD/Const::QE);
    // calculate the un-normalized electric potential energy
    pe = 0.5*Const::EPS_0*(eft*eft);
    
    // Calculate the total thermal energy of all the electrons
    double Th = (tempE*Const::eV)*(species.spwt)*species.numparticle;
    // normalize the potential energy by the total thermal 
    pe = pe/Th;
    return pe;
}

