#ifndef _INIT_H_
#define _INIT_H_

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "species.h"
#include "domain.h"


class Init 
{
public:
    //in the constructor we pass Species and Domain class instances to access the species and domain class
    //memebr variable and methods 
    Init(Species &species, Domain &domain);
    void display();

private:
    Species &species;
    Domain &domain;

    double SampleVel(Species &species);


};

#endif 
