#include "species.h"


/* note :
"-> " symbol means dereferencing a pointer , here we created pointers to instances of class named domain.
so to dereference the instances we have to use the symbol "->"*/

Species::Species(string name, double mass, double charge, double spwt, double temp, int numparticle, double vs,double fract_den , std:: string initialization, Domain &domain):domain(domain)
{
    this->name = name;
    this->mass = mass;
    this->charge = charge;
    this->spwt = spwt;
    this->temp = temp;
    this->numparticle = numparticle;
    this->vs = vs;
    this->fract_den = fract_den;
    this->initialization = initialization;

    den = vec<double>(domain.ni);

    velmesh = vec<double>(domain.ni);
    
    //velprev = vec<double>(numparticle);

    for(int i = 0; i < domain.num_threads; i++)
    {
        domain.buffers.emplace_back(domain.ni);
    }
    
}

void Species::AddParticle(Particle part)
{
    part_list.push_back(part);
}

void Species::Push_species_serial(vector<Species> &species_list, int sub_cycle)
{

    for(auto it = part_list.begin(); it != part_list.end();)
    {
        Particle &part = *it;
        //double lc = domain.XtoL(part.pos);
        double lc = part.pos;
        //cout<<part.pos<<endl;
        // compute charge to mass ratio
        double qm = charge/mass;

        double part_ef = domain.Interpolate(lc,domain.ef);
        //cout<<part_ef<<endl;
        double wl = domain.LDe*domain.LDe*domain.wpe*domain.wpe;

        part.vel[0] += qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_ef*(domain.DT*sub_cycle);

        //Advance particle position

        double old_pos = part.pos;

        //part.pos += ((domain.wp*domain.LD)/(domain.L*domain.W))*part.vel*(domain.DT*sub_cycle);
        part.pos += ((domain.vel_norm)/(domain.L*domain.W))*part.vel[0]*(domain.DT*sub_cycle);

        double new_pos = part.pos;

        if(fabs(new_pos - old_pos) >= 1)
        {
            if(name  == "electron")
            {
                domain.ele_cross++;
                domain.crossed_cellno_ele += int(fabs(new_pos - old_pos));
            }
            if(name  == "ion")
            {
                domain.ion_cross++;
                domain.crossed_cellno_ion += int(fabs(new_pos - old_pos));
            }
        }

        if(domain.bc == "pbc")
        {
            if (part.pos < domain.x0)
		    {
		        part.pos = part.pos + domain.xL;
		    }
		    else if(part.pos >= domain.x0 + domain.xL)
		    {
			    part.pos = part.pos - domain.xL;
		    }

            //display::print(part.pos);
        }

        else if (domain.bc == "rbc")
        {
            if (part.pos < domain.x0)
            {
                part.pos = domain.x0 + (domain.x0 - part.pos); // Reflect position
                part.vel[0] = -part.vel[0]; // Reverse velocity
            }
            else if (part.pos >= domain.x0 + domain.xL)
            {
                part.pos = domain.x0 + domain.xL - (part.pos - (domain.x0 + domain.xL)); // Reflect position
                part.vel[0] = -part.vel[0]; // Reverse velocity
            }
        }

        else if(domain.bc == "open")
        {
    
            if (part.pos < domain.x0 || part.pos >= domain.x0 + domain.xL)
            {
                if(part.pos < domain.x0)
                {
                    domain.wall_left = true;
                    if(IsIon())
                    {
                        SEE(species_list[0],domain);
                    }

                    // Accumulate charge at the left wall
                    if (IsIon())
                    {
                        domain.vL += 1*spwt/domain.density;  // Add ion charge
                    }
                    else
                    {
                        domain.vL += -1*spwt/domain.density ; 
                    }
                }
                if(part.pos >= domain.x0 + domain.xL)
                {
                    domain.wall_left = false;
                    if(IsIon())
                    {
                        SEE(species_list[0],domain);
                    }
                    //new code line
                    // Accumulate charge at the right wall
                    if (IsIon())
                    {
                        domain.vR += 1*spwt/domain.density;  // Add ion charge
                    }
                    else
                    {
                        domain.vR += -1*spwt/domain.density; 
                    }
                }

                
                // Compute wall currents after charge accumulation
                double J_L = domain.vL / (domain.DT * sub_cycle);
                double J_R = domain.vR / (domain.DT * sub_cycle);
                
                domain.I_leftwall = J_L;
                domain.I_rightwall = J_R;

                it = part_list.erase(it);
                continue;
            }
        }
        ++it;

    }
}

void Species::update(int start, int end, int sub_cycle, vector<Species> &species_list)
{
    //Use iterators for parallel processing
    auto it = part_list.begin() + start;
    auto end_it = part_list.begin() + end;
    
    int index = 0;

    while (it != end_it)
    {
        Particle &part = *it;
        // double lc = domain.XtoL(part.pos);
        double lc = part.pos;
        // Compute charge to mass ratio
        double qm = charge / mass;

        // Interpolate electric field at particle location
        double part_ef = domain.Interpolate(lc, domain.ef);
        double wl = domain.LDe * domain.LDe * domain.wpe * domain.wpe;

        //velprev(index++) = part.vel;
        part.prevvel[0] = part.vel[0];

        part.prevvel[1] = part.vel[1];
        part.prevvel[2] = part.vel[2];
        
        // Update particle velocity
        part.vel[0] += qm * ((domain.density * Const::QE * domain.L) / (Const::EPS_0 * domain.W * domain.vel_norm)) * part_ef * (domain.DT * sub_cycle);

        // Advance particle position
        double old_pos = part.pos;
        part.pos += ((domain.vel_norm) / (domain.L * domain.W)) * part.vel[0] * (domain.DT * sub_cycle);
        double new_pos = part.pos;
    
        // Check for particle crossing
        if (fabs(new_pos - old_pos) >= 1)
        {
            if (name == "electron")
            {
                domain.ele_cross++;
                domain.crossed_cellno_ele += int(fabs(new_pos - old_pos));
            }
            if (name == "ion")
            {
                domain.ion_cross++;
                domain.crossed_cellno_ion += int(fabs(new_pos - old_pos));
            }
            // cerr << "Warning: Particle crossing a full cell!" << endl;
        }

        // Handle boundary conditions
        if (domain.bc == "pbc")
        {
            if (part.pos < domain.x0)
            {
                part.pos += domain.xL;
            }
            else if (part.pos >= domain.x0 + domain.xL)
            {
                part.pos -= domain.xL;
            }
        }

        else if (domain.bc == "rbc")
        {
            if (part.pos < domain.x0)
            {
                part.pos = domain.x0 + (domain.x0 - part.pos); // Reflect position
                part.vel[0] = -part.vel[0]; // Reverse velocity
            }
            else if (part.pos >= domain.x0 + domain.xL)
            {
                part.pos = domain.x0 + domain.xL - (part.pos - (domain.x0 + domain.xL)); // Reflect position
                part.vel[0] = -part.vel[0]; // Reverse velocity
            }
        }

        else if(domain.bc == "open")
        {
            if (part.pos < domain.x0 || part.pos >= domain.x0 + domain.xL)
            {
                if(part.pos < domain.x0)
                {
                    domain.wall_left = true;
                    if(IsIon())
                    {
                        SEE(species_list[0],domain);
                    }

                    // Accumulate charge at the left wall
                    if (IsIon())
                    {
                        domain.vL += 1*spwt/domain.density;  // Add ion charge
                    }
                    else
                    {
                        domain.vL += -1*spwt/domain.density ; 
                    }
                }
                if(part.pos >= domain.x0 + domain.xL)
                {
                    domain.wall_left = false;
                    if(IsIon())
                    {
                        SEE(species_list[0],domain);
                    }
                    //new code line
                    // Accumulate charge at the right wall
                    if (IsIon())
                    {
                        domain.vR += 1*spwt/domain.density;  // Add ion charge
                    }
                    else
                    {
                        domain.vR += -1*spwt/domain.density; 
                    }
                }
                it = part_list.erase(it);
                continue;
            }
        }
        ++it;
    }
}

void Species::Push_species(vector<Species> &species_list, int sub_cycle)
{
    
    if(domain.push_parallal)
    {
        int particles_per_threads = numparticle/domain.num_threads;
       
        std::vector<thread> threads;

        for (int i = 0; i < domain.num_threads; ++i) 
        {
            int start = i * particles_per_threads;
            int end = (i == domain.num_threads - 1) ? part_list.size() : (i + 1) * particles_per_threads;

            threads.emplace_back(&Species::update, this, start, end, sub_cycle, std::ref(species_list));
        }

        // Join threads to wait for their completion
        for (auto &thread : threads) 
        {
            thread.join();
        }
    }
    else
    {
        Push_species_serial( species_list, sub_cycle);
    }
}

void Species::ScatterSpecies_serial()
{   
    //reset density
    den = 0;
    for(Particle &part : part_list)
    {
        double lc = domain.XtoL(part.pos);
        domain.Scatter(lc,spwt,den);
    }

    den /= (domain.dx*domain.L);
    den /= domain.density;

    if(domain.bc == "pbc")
    {
        den(0) += den(domain.ni-1);
        den(domain.ni-1) =  den(0);
    }
    else if(domain.bc == "open" || domain.bc == "rbc")
    {
        den(0) *= 2;
        den(domain.ni-1) *= 2;
    }
}

/// mesh averaged velocity
void Species::ScatterVel_serial()
{   
    //reset density
    velmesh = 0;
    for(Particle &part : part_list)
    {
        double lc = domain.XtoL(part.pos);
        domain.Scatter(lc,spwt*part.vel[0],velmesh);
    }

    velmesh /= (domain.dx*domain.L);

    for(int i = 0 ; i < domain.ni ; i++)
    {
        velmesh(i) /= den(i)*domain.density;
    }
    //velmesh /= den*domain.density;

    if(domain.bc == "pbc")
    {
        velmesh(0) += velmesh(domain.ni-1);
        velmesh(domain.ni-1) =  velmesh(0);
    }
    else if(domain.bc == "open" || domain.bc == "rbc")
    {
        velmesh(0) *= 2;
        velmesh(domain.ni-1) *= 2;
    }
}

void Species::parall_deposit(int threadidx, int start, int end)
{
    domain.buffers[threadidx].clear();
    for(int i = start; i < end; i++)
    {
        Particle &part = part_list[i];
        double lc = domain.XtoL(part.pos);
        domain.Scatter(lc,spwt,domain.buffers[threadidx]);
    }
}


void Species::ScatterSpecies()
{
    if(domain.deposit_parallal)
    {
        int particles_per_threads = numparticle/domain.num_threads;
       
        std::vector<thread> threads;

        domain.buffers.resize(domain.num_threads);
        
        for (int i = 0; i < domain.num_threads; i++) 
        {
            int start = i * particles_per_threads;
            int end = (i == domain.num_threads - 1) ? part_list.size() : start + particles_per_threads;
            threads.emplace_back(&Species::parall_deposit,this,i, start, end);
        }

        //Join threads 
        for (auto &thread : threads) 
        {
            thread.join();
        }
        
        //reset density
        den.clear();

        for(int i = 0 ; i < domain.num_threads; i++)
        {
            den += domain.buffers[i];
        }

        den /= (domain.dx*domain.L);
        den /= domain.density;

        if(domain.bc == "pbc")
        {
            den(0) += den(domain.ni-1);
            den(domain.ni-1) =  den(0);
        }
        else if(domain.bc == "open" || domain.bc == "rbc")
        {
            den(0) *= 2;
            den(domain.ni-1) *= 2;
        }
    }
    else
    {
        ScatterSpecies_serial();
    }
}

void Species::Rewind_species()
{
    for (Particle &part: part_list)
    {
        double lc = domain.XtoL(part.pos);
        //double lc = part.pos;
        // compute charge to mass ratio
        double qm = charge/mass;
        //cout<<qm<<endl;
        double part_ef = domain.Interpolate(lc,domain.ef);
        double wl = domain.LDe*domain.LDe*domain.wpe*domain.wpe;

        part.vel[0] -= 0.5*qm*((domain.density*Const::QE*domain.L)/(Const::EPS_0*domain.W*domain.vel_norm))*part_ef*domain.DT;
                
    }
}

double Species::Compute_KE(Species &species)
{
    double ke = 0.0;

    for (Particle &part : part_list)
    {
        // Compute average velocity in each direction
        double v_avg_x = (part.vel[0]);
        double v_avg_y = (part.vel[1]);
        double v_avg_z = (part.vel[2]);

        // Compute squared velocity magnitude
        double v_squre = v_avg_x * v_avg_x + v_avg_y * v_avg_y + v_avg_z * v_avg_z;

        // Kinetic energy: KE = 0.5 * m * v^2
        ke += v_squre * (domain.vel_norm * domain.vel_norm);
    }

    // Multiply by 0.5 * mass * spwt (statistical weight)
    ke *= 0.5 * (spwt * mass);

    // Determine normalization factor based on normscheme
    double Th;
    if (domain.normscheme == 5)
    {
        Th = domain.energy_scale;
    }
    else
    {
        Th = (species.temp * Const::eV) * (species.spwt) * species.numparticle;
    }

    // Normalize kinetic energy
    ke /= Th;
    return ke;
}

std::tuple<double, double, double> Species::Compute_Momentum(Species &species)
{
    double p_x = 0.0, p_y = 0.0, p_z = 0.0;

    for (Particle &part : part_list)
    {
        // Compute total momentum in each direction (average velocity)
        p_x += (part.vel[0]) * domain.vel_norm;
        p_y += (part.vel[1])* domain.vel_norm;
        p_z += (part.vel[2])* domain.vel_norm;
    }

    // Multiply by mass and spwt (statistical weight)
    p_x *= (spwt * mass);
    p_y *= (spwt * mass);
    p_z *= (spwt * mass);

    // Normalize momentum
    double Thp;
    if (domain.normscheme == 5)
    {
        Thp = sqrt(domain.energy_scale * species.spwt * species.numparticle * Const::ME);
    }
    else
    {
        Thp = sqrt((species.temp * Const::eV) / Const::ME) * species.spwt * species.numparticle * Const::ME;
    }

    p_x /= Thp;
    p_y /= Thp;
    p_z /= Thp;

    return std::make_tuple(p_x, p_y, p_z);
}

void SEE(Species &species, Domain &domain)
{
    int num_see = int(domain.see_rate + rnd());
    if(domain.see_rate == 0)
    {
        num_see = 0;
    }
    double x,v;
    for(int i = 0; i < num_see; i++)
    {
        if(domain.wall_left)
        {
            x = domain.x0;
            v = Init::SampleVel(species,domain.tempwall);
        }
        else
        {
            x = domain.xL-rnd();
            v = - Init::SampleVel(species,domain.tempwall);
        }

        species.AddParticle(Particle(x,v,0,0));
    }
}

bool Species::IsIon()
{
    if(charge > 0)
    {return true;}
    else
    {return false;}
}