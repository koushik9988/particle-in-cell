#include "species.h"
#include <vector>
#include "domain.h"
#include <cstring>

/* note :
"-> " symbol means dereferencing a pointer , here we created pointers to instances of class named domain.
so to dereferenced the instances we have used the symbol "->"*/

Species::Species(string name, double mass, double charge, double spwt, double temp, int numparticle,Domain &domain):domain(domain)
{
    this->name = name;
    this->mass = mass;
    this->charge = charge;
    this->spwt = spwt;
    this->temp = temp;
    this->numparticle = numparticle;

    //Domain *domain;
    den =  new double[domain.ni];
    memset(den, 0, domain.ni*sizeof(double));
}

/*
Species::~Species() 
{   
    std::cout<<"species destructor called"<<name<<std::endl;
    // Deallocate memory for den
    delete[] den;
}*/

void Species::AddParticle(Particle part)
{
    part_list.push_back(part);
}

void Species::Push_species()
{
    for (Particle &part: part_list)
    {
        double lc = domain.XtoL(part.pos);
        //cout<<part.pos<<endl;
        // compute charge to mass ratio
        double qm = charge/mass;

        double part_ef = domain.Interpolate(lc,domain.ef);
        //cout<<part_ef<<endl;
        double wl = domain.LD*domain.LD*domain.wp*domain.wp;

        part.vel += (1/wl)*((qm*domain.tempE*Const::eV)/Const::QE)*part_ef*domain.DT;
        //cout<<domain.DT<<endl;
        //part.vel += qm*(massE/e)*part_ef*DT;
        // Advance particle position
        part.pos += part.vel*domain.DT;
        //cout <<"pos:"<<part.pos<<"\t"<<"vel:"<<part.vel<<endl;
        
		if (part.pos < domain.x0)
		{
		    part.pos = part.pos + domain.xL;
		}
		else if(part.pos >= domain.x0 + domain.xL)
		{
			part.pos = part.pos - domain.xL;
		}
    }
}

void Species::update(int start, int end)
{
    for(int i = start; i < end; i++)
    {
        Particle &part = part_list[i];
        double lc = domain.XtoL(part.pos);
        //cout<<part.pos<<endl;
        // compute charge to mass ratio
        double qm = charge/mass;

        double part_ef = domain.Interpolate(lc,domain.ef);
        //cout<<part_ef<<endl;
        double wl = domain.LD*domain.LD*domain.wp*domain.wp;

        part.vel += (1/wl)*((qm*domain.tempE*Const::eV)/Const::QE)*part_ef*domain.DT;
        //cout<<domain.DT<<endl;
        //part.vel += qm*(massE/e)*part_ef*DT;
        // Advance particle position
        part.pos += part.vel*domain.DT;
        //cout <<"pos:"<<part.pos<<"\t"<<"vel:"<<part.vel<<endl;
        
		if (part.pos < domain.x0)
		{
		    part.pos = part.pos + domain.xL;
		}
		else if(part.pos >= domain.x0 + domain.xL)
		{
			part.pos = part.pos - domain.xL;
		}

    }
}

void Species::Push_species_parallel()
{
    int num_threads = std::thread::hardware_concurrency();
    int particles_per_threads = numparticle/num_threads;
    //std::cout<<particles_per_threads<<std::endl;
    //int particles_per_threads = part_list.size()/num_threads;

    std::vector<thread> threads;

    for (int i = 0; i < num_threads; ++i) 
    {
        int start = i * particles_per_threads;
        int end = (i == num_threads - 1) ? part_list.size() : (i + 1) * particles_per_threads;

        threads.emplace_back(&Species::update, this, start, end);
    }

    // Join threads to wait for their completion
    for (auto &thread : threads) 
    {
        thread.join();
    }

}

void Species::ScatterSpecies()
{   
    double *field = den;
    //memset(field,0,sizeof(double)*domain.ni);

    for(Particle &part : part_list)
    {
        double lc = domain.XtoL(part.pos);
        //std::cout << "Particle Position: " << part.pos << ", lc: " << lc << std::endl;
        //std::cout<<lc<<std::endl;
        domain.Scatter(lc,spwt,field);
    }

    for(int i=0; i<domain.ni; i++)
    {
    	field[i] /= domain.dx;
    }

	/*Normalize the field value*/
	for(int i=0; i<domain.ni; i++)
    {
		field[i] /= domain.density;
        //cout<<den[i]-field[i]<<endl;
    }

    field[0] += field[domain.ni-1];
    field[domain.ni-1] =  field[0];
    
    //domain.printmat(field);
}


void Species::Rewind_species()
{
    for (Particle &part: part_list)
    {
        double lc = domain.XtoL(part.pos);
        // compute charge to mass ratio
        double qm = charge/mass;
        //cout<<qm<<endl;
        double part_ef = domain.Interpolate(lc,domain.ef);
        double wl = domain.LD*domain.LD*domain.wp*domain.wp;

        //part.vel -= 0.5*(1/wl)*Const::QE*(qm*domain.tempE/Const::QE)*part_ef*domain.DT;
        part.vel -= 0.5*(1/wl)*((qm*domain.tempE*Const::eV)/Const::QE)*part_ef*domain.DT;
        //cout<<part.vel<<endl;
        
    }
}

double Species::Compute_KE(Species &species)
{
    double ke = 0;
    for (Particle &part:part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke += (part.vel*part.vel)*(domain.wp*domain.LD)*(domain.wp*domain.LD);
    }
    /*Multiply 0.5*mass for all particles*/
    ke *= 0.5*(spwt*mass);
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (species.temp*Const::eV)*(species.spwt)*species.numparticle;

    // Normalize the kinetic energy by the total thermal energy of cold electrons    
    ke = ke/Th;
    return ke;
}