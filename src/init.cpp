#include "init.h"

Init::Init(Species &species, Domain &domain) : species(species),domain(domain)
{

    for(int p = 0; p < species.numparticle ; p++) 
    {
        double x = 0;
        double v = 0;

        double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);

        if(species.name == "electron") 
        {
            x = domain.x0 + p * ((domain.xL) / species.numparticle);
            //double v_th = sqrt(2 * Const::K_b * species.temp * Const::EV_to_K / species.mass);
            //x = domain.x0 + (domain.xL)*rnd();
            v = SampleVel(species) + domain.v_e*v_th;
            v = v/v_th;      
        } 
        else if(species.name == "ion") 
        {
            x = domain.x0 + p * ((domain.xL) / species.numparticle);
            //double v_th = sqrt(2 * Const::K_b * species.temp * Const::EV_to_K / species.mass);
            //x = domain.x0 + (domain.xL)*rnd();
            v = SampleVel(species) + domain.v_i*v_th;
            v = v/v_th;
        }

        else if(species.name == "negion") 
        {
            x = domain.x0 + p * ((domain.xL) / species.numparticle);
            //double v_th = sqrt(2 * Const::K_b * species.temp * Const::EV_to_K / species.mass);
            //x = domain.x0 + (domain.xL)*rnd();
            v = SampleVel(species) + domain.v_n*v_th;
            v = v/v_th;
        }

        else if(species.name == "beam") 
        {
            x = domain.x0 + p * ((domain.xL) / species.numparticle);
            //x = 3*domain.xL/4;
            double v_thb = sqrt(2 * Const::K_b * species.temp * Const::EV_to_K / species.mass);
            v = SampleVel(species) + domain.v_b*v_thb;
            v = v/v_th;
        }

        // Add to the particle list/vector.
        species.AddParticle(Particle(x, v));
    }
}

double Init::SampleVel(Species &species)
{
    //double v_th = sqrt(2 * Const::K_b * domain.tempE * Const::EV_to_K / Const::ME);
    double v_th =  sqrt(2*Const::K_b*species.temp*Const::EV_to_K/species.mass);
    double vt = v_th * sqrt(2) * (rnd() + rnd() + rnd() - 1.5) ;//+ domain.v_i*domain.wp * domain.LD;
    return vt;

}
void Init::display()
{
    cout<<species.numparticle<<endl;
}
