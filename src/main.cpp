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

int main() 
{     
    //parsing input.ini file and storing values
    const std::string filename = "input.ini";
    auto iniData = INIParser::parse(filename);

    std::string outputfolder = INIParser::getString(iniData["file"],"output");
    //grid or domain
    int NC = INIParser::getInt(iniData["domain"], "NC");//cell no
    int ni = NC+1; // no of grid points is one more than cell no
    //simulation
    double v_i = INIParser::getDouble(iniData["simulation"], "v_i");
    double v_e = INIParser::getDouble(iniData["simulation"], "v_e");
    double v_b = INIParser::getDouble(iniData["simulation"], "v_b");
    double v_n = INIParser::getDouble(iniData["simulation"], "v_n");
    int NUM_TS = INIParser::getInt(iniData["time"], "NUM_TS");
    double DT_coeff = INIParser::getDouble(iniData["diagnostics"],"DT_coeff");
    int save_fig = INIParser::getDouble(iniData["diagnostics"],"save_fig");
    int write_interval = INIParser::getInt(iniData["diagnostics"],"write_interval");
    int write_interval_phase = INIParser::getInt(iniData["diagnostics"],"write_interval_phase");
    int write_diagnostics = INIParser::getInt(iniData["diagnostics"],"write_diagnostics");
    double x0 = INIParser::getDouble(iniData["domain"],"x0");
    double alpha = INIParser::getDouble(iniData["simulation"],"alpha");
    double beta = INIParser::getDouble(iniData["simulation"],"beta");
    int write_flag = INIParser::getInt(iniData["diagnostics"],"write_flag");
    //population
    int nE = INIParser::getInt(iniData["population"], "nParticlesE");
    int nI = INIParser::getInt(iniData["population"], "nParticlesI");
    int nB = INIParser::getInt(iniData["population"], "nParticlesB");
    int nN = INIParser::getInt(iniData["population"], "nParticlesN");
    double tempE = INIParser::getDouble(iniData["population"], "tempE");
    double tempI = INIParser::getDouble(iniData["population"], "tempI");
    double tempN = INIParser::getDouble(iniData["population"], "tempN");
    double tempB = INIParser::getDouble(iniData["population"], "tempB");
    double massE = Const::ME;
    double massI = INIParser::getDouble(iniData["population"], "massI");
    massI = Const::AMU*massI;
    double massB = INIParser::getDouble(iniData["population"], "massB");
    massB = Const::AMU*massB;
    double massN = INIParser::getDouble(iniData["population"], "massN");
    massN = Const::AMU*massN;
    double den = INIParser::getDouble(iniData["simulation"],"density");
    double species_no = INIParser::getInt(iniData["simulation"],"number_of_species");
    std:: string bc = INIParser::getString(iniData["simulation"],"bc");
    
    //normalized density 
    double ni0 = den;
    double ne0 = ni0*(1- alpha - alpha*beta);
    double nn0 = alpha*ni0; //fraction of negative ion to background positive ion(eletronegativity)*/
    double nb0 = beta*nn0; //fraction of negative beam to background negative ion 

    //noramlizing quantity(electron)
    double LD = sqrt((Const::EPS_0*Const::K_b*tempE*Const::EV_to_K)/(ne0*Const::QE*Const::QE)); // Electron Debye Length    
    double wp = sqrt((ne0*Const::QE*Const::QE)/(massE*Const::EPS_0)); // Total Electron Plasma Frequency
    double wpi = sqrt((ni0*Const::QE*Const::QE)/(massI*Const::EPS_0)); //ion timescale
    double CS = sqrt(tempE*Const::K_b*Const::EV_to_K/massI); // Ion acoustic speed

    print("ion-acosutic speed :", CS);
    print("debye lenght :", LD);

    double stepSize =  LD;
    //normalized spacing
    double dx = stepSize/LD;
    
    double DT = DT_coeff*(1.0/wp);
    //normalized time step
    DT = wp*DT;

    print("normalized cell spacing :", dx);
    print("normalized time step :",DT);
    double f = (wp / (2 * Const::PI));
    print("time period :",1.0/f);
    print("ion scale :",1.0/(wpi/(2*Const::PI)));
    //print("time period 2 :",1.0/wp);

    Domain domain(x0,dx,ni);

    domain.bc  = bc;
    domain.set_normparam(LD,wp);
    domain.set_simparam(tempE, tempI, tempN,tempB, den, v_e, v_i, v_n, v_b, alpha,species_no);
    domain.set_time(DT,NUM_TS,write_interval);

    domain.display();

    Output output(outputfolder,domain);

    double electron_spwt = ne0*domain.xL*LD/nE;
    double ion_spwt = ni0*domain.xL*LD/nI;
    double negion_spwt = nn0*domain.xL*LD/nN;
    double beam_spwt = nb0*domain.xL*LD/nB;

    vector <Species> species_list;
    
    species_list.emplace_back("electron", massE,-Const::QE, electron_spwt, tempE, nE, domain);
    species_list.emplace_back("ion", massI, Const::QE, ion_spwt, tempI, nI, domain);
    species_list.emplace_back("negion", massN, -Const::QE, negion_spwt, tempN, nN, domain);
    species_list.emplace_back("beam", massB, -Const::QE, beam_spwt, tempB, nB, domain);


    output.write_metadata(NC,NUM_TS,write_interval,write_interval_phase,DT_coeff, nE, nI, nN, nB,
     tempE,tempI,tempB, tempN, alpha, beta, massI, massN, massB,den,save_fig,v_e, v_i,v_n,v_b);

    for(Species &p:species_list)
    {
         cout<< p.name << '\n' << p.mass<< '\n' << p.charge << '\n' << p.spwt << '\n' << p. numparticle <<'\n'<<p.charge_sig<<endl;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();

    //initializing the species by creating instances if Init class
    for(Species &sp : species_list)
    {
        Init init(sp,domain);
    } 
    
    FieldSolve fieldsolver(domain);
    
    for (Species &sp:species_list)
	{
		sp.ScatterSpecies();
	}

    domain.ComputeRho(species_list);

    fieldsolver.SolvePotDirect();
    fieldsolver.CalculateEfield();

    for (Species &sp:species_list)
	{
        sp.Rewind_species();
    }

    //int index = 0;

    //--------------MAIN LOOP----------------------- 
    for(int ts = 0 ; ts < NUM_TS + 1; ts++)
    {
        for (Species &sp:species_list)
		{
			sp.ScatterSpecies();
		}

        domain.ComputeRho(species_list);
        
        fieldsolver.SolvePotDirect();
        fieldsolver.CalculateEfield();

        for (Species &sp:species_list)
		{
			if(domain.bc == "pbc")
            {
                //sp.Push_species();
                sp.Push_species_parallel();
            }
            else if(domain.bc == "open")
            {
                sp.Push_species_open();
                //sp.Push_species_reflect();
            }
            
		}
    
        double max_phi = domain.phi[0];
        for(int i=0; i<domain.ni; i++)
            if (domain.phi[i]>max_phi) max_phi = domain.phi[i];
                
        if(ts%write_interval== 0)
        {
            //int k = int(index/write_interval);
            if(write_flag == 1 || write_flag == 2 )
            {    
                output.write_field_data(ts);
                output.storeKE_to_matrix(ts,species_list);
                for(Species &sp:species_list)
                {
                    output.write_den_data(ts,sp);
                }    
            }
        }
        
        if(ts%write_interval_phase == 0)
        {
            if(write_flag == 1 || write_flag == 3)
            {
                for(Species &sp:species_list)
                {
                    output.write_particle_data(ts,sp);
                }
            }
        }
        
        if(ts%write_diagnostics == 0)
        {
            output.diagnostics(ts, max_phi,species_list);
        }  
    }

    //output.printmatrix(int(NUM_TS/write_interval) +1 ,species_list.size() + 1,output.store_ke);

    output.write_ke();

    auto end_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    std::cout << "Elapsed time: " << elapsed_time.count()/60 << "minutes." << std::endl;

    return 0;
}
