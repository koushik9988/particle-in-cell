#include "init.h"

Init::Init(Species &species, Domain &domain) : species(species),domain(domain)
{

    auto initilization = INIParser::loadtypeextract(species.initialization);

    auto init_type = initilization.first;
    int n = initilization.second;

    ///
    // Create a random number generator
    std::random_device rd; 
    std::mt19937 gen(rd());

    double mean = (domain.x0 + domain.xL) / 2.0; 
    double sigma = domain.xL/5;

    std::normal_distribution<> gaussian(mean, sigma);
    ///
    
    for(int p = 0; p < species.numparticle ; p++) 
    {
        double x = 0;
        double vx = 0;
        double vy = 0;
        double vz = 0;

        if(init_type == "random") 
        {
            x = domain.x0 + (domain.xL) * rnd();
            if(x >= domain.xL)
            {
                x = x - domain.xL;
            }
            if(x < 0)
            {
                x = x + domain.xL;
            }
        }
        else if (init_type == "uniform") 
        {
            x = domain.x0 + p * (domain.xL / (species.numparticle -1));
            if(x >= domain.xL)
            {
                x = x - domain.xL;
            }
            if(x < 0)
            {
                x = x + domain.xL;
            }
            //display::print(x);
        }
        else if (init_type == "sin") 
        {
            x = domain.x0 + p * (domain.xL / (species.numparticle -1 ));
            double k = 2*Const::PI*n/domain.xL;
            x = x + sin(k*x);
            if(x >= domain.xL)
            {
                x = x - domain.xL;
            }
            if(x < 0)
            {
                x = x + domain.xL;
            }
        }
        else if (init_type == "cos") 
        {
            x = domain.x0 + p * (domain.xL / (species.numparticle-1));
            double k = 2*Const::PI*n/domain.xL;
            x = x + cos(k*x);
            if(x >= domain.xL)
            {
                x = x - domain.xL;
            }
            if(x < 0)
            {
                x = x + domain.xL;
            }
        }
        else if (init_type == "gaussian") 
        {
            // Sample position from Gaussian distribution
            x = gaussian(gen);

            if (x >= domain.xL)
            {
                x = x - domain.xL;
            }
            if (x < 0)
            {
                x = x + domain.xL;
            }
        }

        vx = SampleVel(species) + species.vs * domain.vel_norm;
        vy = 0;
        vz = 0;

        vx = vx / domain.vel_norm;
        vy = vy / domain.vel_norm;
        vz = vz / domain.vel_norm;

        species.AddParticle(Particle(x, vx,vy,vz));
        //display::print(x);
    }
}

double Init::SampleVel(Species &species)
{
    //double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);
    double v_th =  sqrt(2*Const::K_b*species.temp*Const::EV_to_K/species.mass);
    double vt = v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5) ;//+ domain.v_i*domain.wp * domain.LD;
    //double vt = v_th * sqrt(2) * (rnd()*rnd()*rnd() - 1.5) ;
    return vt;
}

double Init::SampleVel(Species &species, double temp)
{
    //double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);
    double v_th =  sqrt(2*Const::K_b*temp*Const::EV_to_K/species.mass);
    double vt = v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5); //+ domain.v_i*domain.wp * domain.LD;
    return vt;
}
