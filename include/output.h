#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <fstream>
#include <vector>
#include <filesystem>
#include <string>
#include <sstream>
#include "domain.h"
#include "species.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <map>

#include "H5Cpp.h"

using namespace H5;

class Domain;
class Species;

class Output 
{
    public:

    std::map<std::string, Group> particle_groups;
    std::map<std::string, Group> den_subgroups;

    Output(const std::filesystem::path& outputfolder, Domain& domain);
    //~Output();
    void write_particle_data(int ts, Species& species);
    //void write_particle_data(H5::Group& group, int ts, Species& species);
    void write_den_data(int ts,  Species& species);
    void write_field_data(int ts);
    //void write_ke(int ts,std::vector<Species> &species_list);
    void write_ke();
    void storeKE_to_matrix(int ts, std::vector<Species> &species_list);
    void printmatrix(int row, int col, double **matrix);
    void diagnostics(int ts, double max_phi, std::vector<Species> &species_list);
    void write_metadata(int NC, int NUM_TS, int write_int, int write_int_phase, double DT, int nE, int nI, int nN, int nB, double Te, double Ti, double Tb, double Tn, double alpha,
 double beta, double mI, double mN, double mB, double density, int save_fig, double ve, double vi, double vn, double vb);
    double **store_ke;
    int sp_no ;//= species_list.size();
    int t_step;// = int(domain.NUM_TS/domain.write_interval) + 1 ;
    private:
    //Species &species;
    std::filesystem::path outputfolder;
    Domain& domain;
    H5File file; // Declare H5::H5File object to handle HDF5 file operations
    Group field_data_group;
    Group time_group;
    Group metadata_group;
    //std::vector<Species> species_list;
   
};

namespace display
{
    template<typename T>
    void print(const T& value) 
    {
        std::cout << value << std::endl;
    }
    // Recursive template function to print multiple arguments
    template<typename T, typename... Args>
    void print(const T& value, Args&&... args) 
    {
        std::cout << value;
        print(std::forward<Args>(args)...); // Recursive call with the remaining arguments
    }

}

#endif
