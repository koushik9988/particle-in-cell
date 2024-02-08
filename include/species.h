#ifndef _SPECIES_H_
#define _SPECIES_H_

#include <iostream>
#include <list>
#include <vector>
#include "domain.h"
#include <thread>

using namespace std;

class Domain;

class Particle
{
    public:
    double pos;
    double vel;

    //particle constructor
    Particle(double x, double v):pos(x), vel(v){};
};

class Species
{
    public:
    string name;
    double mass;
    double charge;
    double spwt;
    double temp;
    double *den;
    int numparticle;
    int charge_sig;

    // the "part_list" is a linked list that holds instances of the Particle class.
    //list/vector is a template class and Particle is template argument(similar to regular int,double etc)
    //specifying type of element it will hold which is instaances of "Particle" class here. 
    //(C++ STL containers such as vector,list can hold generic data type like class etc.)
    //list<Particle> part_list;
    //or
    vector<Particle> part_list;

    //constructor 
    Species(string name, double mass, double charge, double spwt, double temp, int numparticle,Domain &domain);

    //Destructor
    //~Species();

    //declare member functions or methods
    void AddParticle(Particle part);
    void Push_species();
    void update(int start, int end);
    void Push_species_parallel();
    void ScatterSpecies();
    void Rewind_species();
    double Compute_KE(Species &species);

    private:
    Domain &domain;
    
};
#endif 