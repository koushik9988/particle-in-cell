/*
Implementation file for MCC, DSMC collisional models

code problem:
(1) p_coll is very small close to zero ?
(2) Manually setting p_coll between 0 and 1 cause soem nan value problem in velocity.
*/

#include "collision.h"
#include "cross_sections.h"

CollisionHandler::CollisionHandler(Domain &domain) : domain(domain)
{
    sigma_tot_e.resize(CS_RANGES, 0.0);
    sigma_ela.resize(CS_RANGES, 0.0);
    sigma_exc.resize(CS_RANGES, 0.0);
    sigma_ionz.resize(CS_RANGES, 0.0);
}


void CollisionHandler::set_electron_cross_sections()
{
    //std::cout << ">> Setting e- / Ar cross sections" << std::endl;
    std::vector<double> energy(CS_RANGES);
    energy[0] = DE_CS;
    std::generate(energy.begin() + 1, energy.end(), [i = 1]() mutable { return DE_CS * (i++); });

    std::transform(energy.begin(), energy.end(), sigma_ela.begin(),[this](double energy_val) { return compute_elastic_CS(energy_val, domain); });
    
    std::transform(energy.begin(), energy.end(), sigma_exc.begin(),[this](double energy_val) { return compute_excitation_CS(energy_val, domain); });
    
    std::transform(energy.begin(), energy.end(), sigma_ionz.begin(),[this](double energy_val) { return compute_ionization_CS(energy_val, domain); });
}

void CollisionHandler::calc_total_cross_sections()
{
    for (size_t i = 0; i < CS_RANGES; ++i)
    {
        sigma_tot_e[i] = (sigma_ela[i] + sigma_exc[i] + sigma_ionz[i]) * domain.GAS_DENSITY;
    }
}

double CollisionHandler::max_electron_coll_freq()
{
    double e, v, nu, nu_max = 0.0;
    for (int i = 0; i < CS_RANGES; ++i)
    {
        e = i * DE_CS;
        v = sqrt(2.0 * e * Const::eV / Const::ME);
        nu = v * sigma_tot_e[i];
        if (nu > nu_max)
        {
            nu_max = nu;
        }
    }
    return nu_max;
}

void CollisionHandler::collision_electron(double xe, double &vxe, double &vye, double &vze, int eindex, Species &species1, Species &species2)
{

    const double F1 = Const::ME / (Const::ME + species2.mass);
    const double F2 = species2.mass / (Const::ME + species2.mass);
    double  t0, t1, t2;//rnd;
    double  g, g2, gx, gy, gz, wx, wy, wz, theta, phi;
    double  chi, eta, chi2, eta2, sc, cc, se, ce, st, ct, sp, cp, energy, e_sc, e_ej;
    double g_before, g_after;

    // Calculate relative velocity before collision & velocity of the center of mass before collision
    gx = vxe;
    gy = vye;
    gz = vze;
    g  = sqrt(gx * gx + gy * gy + gz * gz);
    g_before = g;
    wx = F1 * vxe;
    wy = F1 * vye;
    wz = F1 * vze;

    //display::print("gbefore:",g); //gx,gy,gz giving nan value have to fix.

    // Find Euler angles
    if (gx == 0) 
    {
        theta = 0.5 * Const::PI;
    }
    else 
    {
        theta = atan2(sqrt(gy * gy + gz * gz),gx);
    }
    if (gy == 0)
    {
        if (gz > 0)
        {
            phi = 0.5 * Const::PI;
        }
        else
        {
            phi = - 0.5 * Const::PI;
        }
    }
    else
    {
        phi = atan2(gz, gy);
    }
    st  = sin(theta);
    ct  = cos(theta);
    sp  = sin(phi);
    cp  = cos(phi);

    // Choose the type of collision based on cross-sections
    t0 = sigma_ela[eindex];
    t1 = t0 + sigma_exc[eindex];
    t2 = t1 + sigma_ionz[eindex];
    //rnd = rand();  // Random number between 0 and 1

    if (rnd() < (t0 / t2))
    {  
        //display::print("elatic!");
        // Elastic scattering
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2*Const::PI * rnd();
        
    }
    else if (rnd() < (t1 / t2))
    {  
        //display::print("inelastic!");
        // Excitation
        energy = 0.5 * Const::ME *g * g * domain.vel_norm * domain.vel_norm;  // Energy in joules
        energy = fabs(energy - E_EXC_TH * Const::eV);  // Energy loss for excitation (exitation enegry converted to joule)
        g = sqrt(2.0 * energy / Const::ME);  // Relative velocity after energy loss
        g = g / domain.vel_norm;  // Normalize the velocity
        chi = acos(1.0 - 2.0 * rnd());
        eta = 2*Const::PI * rnd();
    }
    else
    {  
        // Ionization
        energy = 0.5 * Const::ME * g * g * domain.vel_norm * domain.vel_norm; // Energy in joules
        energy = fabs(energy - E_ION_TH * Const::eV);  // Energy loss for ionization
        e_ej = 10.0 * tan(rnd() * atan(energy / Const::eV / 20.0)) * Const::eV;  // Emitted electron energy
        e_sc = fabs(energy - e_ej);  // Incoming electron energy after collision
        g = sqrt(2.0 * e_sc / Const::ME);
        g2 = sqrt(2.0 * e_ej / Const::ME);
        g = g / domain.vel_norm;  // Normalize the velocity
        g2 = g2 / domain.vel_norm;  // Normalize the velocity
        chi = acos(sqrt(e_sc / energy));  // Scattering angle for incoming electron
        chi2 = acos(sqrt(e_ej / energy));  // Scattering angle for emitted electron
        eta = 2*Const::PI * rnd();
        eta2 = eta + Const::PI;

        // Compute velocity of emitted electron
        sc = sin(chi2);
        cc = cos(chi2);
        se = sin(eta2);
        ce = cos(eta2);
        gx = g2 * (ct * cc - st * sc * ce);
        gy = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
        
        // Add the new electron and ion (normalize this as we have to unnormalize for previous calculations) (velocities are converted to lab gram from com frame here)
        species1.AddParticle(Particle(xe, (wx + F2 * gx), (wy + F2 * gy), (wz + F2 * gz)));
        species2.AddParticle(Particle(xe, Init::SampleVel(species2,species2.temp),  Init::SampleVel(species2,species2.temp),  Init::SampleVel(species2,species2.temp)));
    }

    // Scatter the primary electron
    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);

    // Compute new relative velocity 
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    // Post-collision velocity of the colliding electron in lab frame
    vxe = wx + F2 * gx;
    vye = wy + F2 * gy;
    vze = wz + F2 * gz;

    g_after = sqrt(vxe * vxe + vye * vye + vze * vze);

    domain.delta_g = (fabs(g_after - g_before));  // Change in velocity after collision
    //printf("g_before: %f, g_after: %f, delta_g: %f\n", g_before*domain.vel_norm, g_after*domain.vel_norm, domain.delta_g);

}

void CollisionHandler::handle_collisions(Species &electron, Species &target_gas)
{
    for (Particle &part : electron.part_list)
    {
        double v_sqr = (part.vel[0] * part.vel[0] + part.vel[1] * part.vel[1] + part.vel[2] * part.vel[2]) * domain.vel_norm * domain.vel_norm;
        double velocity = sqrt(v_sqr);
        double energy = (0.5 * Const::ME * v_sqr) / Const::eV;
        //int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES - 1);
        int energy_index = std::min(static_cast<int>(energy / DE_CS + 0.5), static_cast<int>(CS_RANGES) - 1);

        double nu = sigma_tot_e[energy_index] * velocity;
        double p_coll = 1 - exp(-nu * (domain.DT/domain.W));
        if (p_coll > rnd())
        {
            collision_electron(part.pos, part.vel[0], part.vel[1], part.vel[2], energy_index, electron, target_gas);
        }
    }
}
