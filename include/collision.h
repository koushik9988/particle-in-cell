#ifndef _COLLISIONAL_MODEL_H_
#define _COLLISIONAL_MODEL_H_

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "domain.h"
#include "species.h"
#include "init.h"

constexpr int CS_RANGES = 1000000; //CS_RANGES is a compile-time constant
constexpr double DE_CS = 0.005;  
//using cross_section = std::array<double, CS_RANGES>;
using cross_section = std::vector<double>;

class CollisionHandler
{
    private:
    Domain &domain;
    cross_section sigma_tot_e;  
    cross_section sigma_ela;
    cross_section sigma_exc;
    cross_section sigma_ionz;
    
    public:
    CollisionHandler(Domain& domain);
    void set_electron_cross_sections();
    void calc_total_cross_sections();
    void collision_electron(double xe, double &vxe, double &vye, double &vze, int eindex, Species &species1, Species &species2);
    void handle_collisions(Species &electron, Species &target_gas);
    double max_electron_coll_freq();
};

#endif  // _COLLISIONAL_MODEL_H_
