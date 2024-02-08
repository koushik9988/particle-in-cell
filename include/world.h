// header file(.h or .hpp) contain function,class,namespace etc declaration
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "species.h"
#include <iostream>
#include <vector>
#include <random>

class Species;

using namespace std;

class Domain
{
    public:
    double x0; //origin points
    double xL; //max lenght of system
    double dx; //cell spacing
    int cell_no; // no of cell
    int ni; //no of grid points
    //double LD;

    //data structures to hold data
    double *phi;
    double *rho;
    double *ef;

    //simulation parameter
    double LD ;
    double wp;
    double tempE;
    double tempI;
    double density;
    double DT;
    double v_i;
    double v_e;

    //constructor 
    Domain(double x0, double dx, int ni);
    //destructor 
    ~Domain();

    void display();

    void set_normparam(double LD, double wp);
    void set_time(double DT);
    void set_simparam(double tempE, double tempI, double density, double  v_e, double v_i);
    /*here we use address of element in the vector container using address we can make changes in the original element 
    here vector<Species> is like kind of int,float etc data type. species is a kind of list/array which elements is different species 
    such as ion,electron and neutral etc.*/
    void ComputeRho(vector<Species> &species); 
    void ComputeRho1(Species &species1, Species &species2);

    //scatter density to mesh/grid.
    void Scatter(double lc, double value , double *field);
    //void ScatterSpecies(Species *species);
    //get density in the grid/mesh points.
    double Interpolate(double lc, double *field);

    double XtoL(double pos);

    void printmat(double *A); 

};

//define a namespace name Const
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	//vacuum permittivity
	const double QE = 1.602176565e-19;		//electron charge
	const double AMU = 1.660538921e-27;		//atomic mass unit
	const double ME = 9.10938215e-31;		//electron mass
	const double K_b = 1.380648e-23;			//Boltzmann constant
	const double PI = 3.141592653;			//pi
	const double EV_to_K = QE/K_b;				//1eV in K ~ 11604
    const double eV = 1.602176565e-19;		//1 eletron volt
    
}

/*object for sampling random numbers*/
class Rnd {
public:
	//constructor: set initial random seed and distribution limits
	Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd is defined somewhere


#endif