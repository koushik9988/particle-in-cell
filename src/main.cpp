/* */
#include "iniparser.h"
#include <iostream>
#include "field.h"
#include "init.h"
#include "domain.h"
#include "species.h"
#include <cmath>
#include <fstream>
#include <time.h>
#include <chrono>
#include "output.h"

using namespace std; 
using namespace display;

FILE *file;

int main() 
{     
    //clock_t start = clock();
    //----------------
    file = fopen("t_den.txt", "w");
    if (!file) 
    {
        std::cerr << "Error opening data file!" << std::endl;
        return 1;
    }
    //----------------

    //parsing input.ini file and storing values
    const std::string filename = "input.ini";
    auto iniData = INIParser::parse(filename);

    
    /* Define Simulation domain*/
    double den;           // Plasma Density
    double stepSize;      // Cell Spacing
    double DT_coeff;
    int write_interval;
    int write_interval_phase;
    double x0;
    //double xL;
    int NC;              // Total number of cells
    int NUM_TS;          // Total Time steps (default)
    int write_flag;
    //int write_interval;
    
    /*Simulation Normalizing Parameters*/
    double ne0, ni0,nb0,nn0;
    double CS;
    double alpha,beta;

    //population
    int nE;     // Number of simulation electrons
    int nI;
    int nB;
    int nN;
    double T_e;     // Temperature of the electrons in eV units
    double T_i;  	  // Temperature of the ions in eV units
    double T_b;
    double T_n;
    double m_e;
    double m_i;
    double m_n;
    double m_b;

    //simulation parameter
    double v_i;      //ion streaming velocity
    double v_e;     //electron streaming velocity
    //double alpha;			 // The fraction of cold electrons
    double v_n;
    double v_b;
    std::string outputfolder;
    
    double dx;
    int ni;

    outputfolder = INIParser::getString(iniData["file"],"output");
    //grid or domain
    NC = INIParser::getInt(iniData["domain"], "NC");//cell no
    ni = NC+1; // no of grid points is one more than cell no

    //simulation
    v_i = INIParser::getDouble(iniData["simulation"], "v_i");
    v_e = INIParser::getDouble(iniData["simulation"], "v_e");
    v_b = INIParser::getDouble(iniData["simulation"], "v_b");
    v_n = INIParser::getDouble(iniData["simulation"], "v_n");
    NUM_TS = INIParser::getInt(iniData["time"], "NUM_TS");
    DT_coeff = INIParser::getDouble(iniData["diagnostics"],"DT_coeff");
    write_interval = INIParser::getInt(iniData["diagnostics"],"write_interval");
    write_interval_phase = INIParser::getInt(iniData["diagnostics"],"write_interval_phase");
    x0 = INIParser::getDouble(iniData["domain"],"x0");
    alpha = INIParser::getDouble(iniData["simulation"],"alpha");
    beta = INIParser::getDouble(iniData["simulation"],"beta");
    write_flag = INIParser::getInt(iniData["diagnostics"],"write_flag");
    //xL

    //Domain domain(x0,dx,ni);
    //population
    nE = INIParser::getInt(iniData["population"], "nParticlesE");
    nI = INIParser::getInt(iniData["population"], "nParticlesI");
    nB = INIParser::getInt(iniData["population"], "nParticlesB");
    nN = INIParser::getInt(iniData["population"], "nParticlesN");
    T_e = INIParser::getDouble(iniData["population"], "tempE");
    T_i = INIParser::getDouble(iniData["population"], "tempI");
    T_n = INIParser::getDouble(iniData["population"], "tempN");
    T_b = INIParser::getDouble(iniData["population"], "tempB");
    m_e = Const::ME;
    m_i = INIParser::getDouble(iniData["population"], "massI");
    m_i = Const::AMU*m_i;
    m_b = INIParser::getDouble(iniData["population"], "massB");
    m_b = Const::AMU*m_b;
    m_n = INIParser::getDouble(iniData["population"], "massN");
    m_n = Const::AMU*m_n;
    den = INIParser::getDouble(iniData["simulation"],"density");
    
    /*
    print("the value of density:",den);
    print("the value of mass of eletron:",m_e);
    print("the value of mass of ion:",m_i);
    print("the value electron temp in ev:",T_e);
    print("the value ion temp in ev:",T_i);
    print("No of eletron:",nE);
    print("No of ion:",nI);*/
    
    /*ni0 = den;
    ne0 = ni0*(1- alpha - beta);
    nb0 = alpha*ni0; //fraction of negative beam to background positive ion
    nn0 = beta*ni0; //fraction of negative ion to background positive ion*/

    
    ni0 = den;
    ne0 = ni0*(1- alpha - alpha*beta);
    nn0 = alpha*ni0; //fraction of negative ion to background positive ion(eletronegativity)*/
    nb0 = beta*nn0; //fraction of negative beam to background negative ion 

    print("electron density",ne0);

    //noramlizing quantity(electron)
    double Debye_lenght = sqrt((Const::EPS_0*Const::K_b*T_e*Const::EV_to_K)/(ne0*Const::QE*Const::QE)); // Electron Debye Length    
    double Plasma_freq = sqrt((ne0*Const::QE*Const::QE)/(m_e*Const::EPS_0)); // Total Electron Plasma Frequency

    //ion normalizing quantity
    //double Debye_lenght = sqrt((Const::EPS_0*Const::K_b*T_i*Const::EV_to_K)/(ni0*Const::QE*Const::QE)); // Electron Debye Length    
    //double Plasma_freq = sqrt((ni0*Const::QE*Const::QE)/(m_i*Const::EPS_0)); // Total Electron Plasma Frequency

    CS = sqrt(T_e*Const::K_b*Const::EV_to_K/m_i); // Ion acoustic speed

    print("ion-acosutic speed :", CS);
    print("debye lenght :", Debye_lenght);

    stepSize = Debye_lenght;
    //normalized spacing
    dx = stepSize/Debye_lenght;
    
    double time_step = DT_coeff*(1.0/Plasma_freq);
    //normalized time step
    time_step = Plasma_freq*time_step;

    print("normalized cell spacing :", dx);
    print("normalized time step :",time_step);
    double f = (Plasma_freq / (2 * Const::PI));
    print("time period :",1.0/f);
    //print(rnd());

    Domain domain(x0,dx,ni);
    
    domain.density = den;
    domain.LD = Debye_lenght;
    domain.wp = Plasma_freq;
    domain.DT = time_step;
    domain.tempE = T_e;
    domain.tempI = T_i;
    domain.tempB = T_b;
    domain.tempN = T_n;
    //domain.alpha = alpha;

    //redefining
    double density = domain.density;
    double LD = domain.LD;
    double wp = domain.wp;
    double DT = domain.DT;
    double tempE = domain.tempE;
    double tempI = domain.tempI;
    double tempN = domain.tempN;
    double tempB = domain.tempB;
    //double alpha = domain.alpha;

    double massE = m_e;
    double massI = m_i;
    double massN = m_n;
    double massB = m_b;
    double nParticlesI = nI;
    double nParticlesE = nE;
    double nParticlesN = nN;
    double nParticlesB = nB;

    domain.set_normparam(LD,wp);
    domain.set_simparam(tempE, tempI, tempN,tempB, density, v_e, v_i, v_n, v_b, alpha);
    domain.set_time(DT);

    domain.display();
    Output output(outputfolder,domain);

    double electron_spwt = ne0*domain.xL/nParticlesE;
    double ion_spwt = ni0*domain.xL/nParticlesI;
    double negion_spwt = nn0*domain.xL/nParticlesN;
    double beam_spwt = nb0*domain.xL/nParticlesB;

    //double TE_e = (tempE*Const::eV)*(electron_spwt)*nParticlesE;

    vector <Species> species_list;
    
    species_list.emplace_back("electron", massE,-Const::QE, electron_spwt, tempE, nParticlesE, domain);
    species_list.emplace_back("ion", massI, Const::QE, ion_spwt, tempI, nParticlesI, domain);
    species_list.emplace_back("negion", massN, -Const::QE, negion_spwt, tempN, nParticlesN, domain);
    species_list.emplace_back("beam", massB, -Const::QE, beam_spwt, tempB, nParticlesB, domain);

    
    
    for(Species &p:species_list)
    {
         cout<< p.name << '\n' << p.mass<< '\n' << p.charge << '\n' << p.spwt << '\n' << p. numparticle <<'\n'<<p.charge_sig<<endl;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();

    //Species &electron = species_list[0];
    //passing address of instances of Species class/ can also use pointers 
    Init init_e(species_list[0],domain), init_i(species_list[1], domain),init_n(species_list[2],domain), init_b(species_list[3], domain);

    FieldSolve fieldsolver(domain);
    
    for (Species &sp:species_list)
	{
		sp.ScatterSpecies();
	}

    domain.ComputeRho(species_list);
    //domain.ComputeRho1(species_list[0],species_list[1]);

    fieldsolver.SolvePotDirect();
    //fieldsolver.SolvePotIter();
    fieldsolver.CalculateEfield();

    for (Species &sp:species_list)
	{
        sp.Rewind_species();
    }
   
    //--------------MAIN LOOP----------------------- 
    for(int ts = 0 ; ts < NUM_TS + 1; ts++)
    {
    
        for (Species &sp:species_list)
		{
			sp.ScatterSpecies();
		}

        domain.ComputeRho(species_list);
        //domain.ComputeRho1(species_list[0],species_list[1]);

        fieldsolver.SolvePotDirect();
        //fieldsolver.SolvePotIter();
        fieldsolver.CalculateEfield();

        for (Species &sp:species_list)
		{
			//sp.Push_species();
            sp.Push_species_parallel();
            //cout<<sp.numparticle<<endl;
		}
        double ke_e = species_list[0].Compute_KE(species_list[0]);
        double ke_i = species_list[1].Compute_KE(species_list[0]);
        double ke_n = species_list[2].Compute_KE(species_list[0]);
        double ke_b = species_list[3].Compute_KE(species_list[0]);
        double pe = domain.ComputePE(species_list[0]);

        //std::cout<<ts<<"\t"<<species_list[0].numparticle<<"\t"<<domain.phi[1]<<endl;
        double max_phi = domain.phi[0];
        for(int i=0; i<domain.ni; i++)
            if (domain.phi[i]>max_phi) max_phi = domain.phi[i];
                
        /*print diagnostics to the screen*/
	    
        fprintf(file, "%f\t%f\n", (ts * DT), species_list[0].den[0]);
        //fprintf(file, "%f\t%f\n", ts * DT / wp, species_list[0].den[20]);
        if(ts%write_interval== 0)
        {
            if(write_flag == 1 || write_flag == 2 )
            {
                output.write_data(ts,species_list);
                output.write_ke(ts,species_list);  
            }
    
        }
        
        if(ts%write_interval_phase == 0)
        {
            printf("TS: %i \t delta_phi: %.3g \t\t nE:%i \t nI:%i \t KE_e:%.3g  \t KE_i:%.3g  \t KE_b:%.3g  \t KE_n:%.3g  \t pe:%.3g\n",
				ts, max_phi-domain.phi[0],species_list[0].numparticle,species_list[1].numparticle,ke_e,ke_i,ke_b,ke_n, pe);

            //output.create_files(ts);
            //output.write_species_data(ts, species_list);
            if(write_flag == 1 || write_flag == 3)
            {
                output.write_particle(ts,species_list[0]);
                output.write_particle(ts,species_list[1]);
                output.write_particle(ts,species_list[2]);
                output.write_particle(ts,species_list[3]);
            }
        }    
    }

    //for(int i =0; i < domain.ni; i++)
    //{
        //fprintf(file, "%f\t%f\n",i*dx, domain.phi[i]);

    //}

    fclose(file);

    //clock_t end = clock();
    auto end_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    std::cout << "Elapsed time: " << elapsed_time.count() << " seconds." << std::endl;

    //cout << "Simulation Time = " << ((end-start)/(double)CLOCKS_PER_SEC)/60 << " minutes" << endl;
    return 0;
}
