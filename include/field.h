#ifndef _FIELD_H_
#define _FIELD_H_

#include "domain.h"
#include <cstring>

class Domain;

class FieldSolve
{
    public:
    //constructor 
    FieldSolve(Domain &domain):domain(domain){};
    
    //void SolvePotDirect(double *x, double *rho);
    void SolvePotDirect();
    bool SolvePotIter();
    //void solvepotiterative();
    //void solvepotcg(); //in future
    //void solvepotspectral();//in future
    //void CalculateEfield(double *phi, double *ef);
    void CalculateEfield();

    protected:
    Domain &domain;

};


#endif