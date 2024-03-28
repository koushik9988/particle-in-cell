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

class Domain;
class Species;

class Output 
{
    public:
    Output(const std::filesystem::path& outputfolder, Domain& domain);
    void write_particle(int ts, Species& species);
    void write_data(int ts, std::vector<Species> &species_list);
    void write_ke(int ts,std::vector<Species> &species_list);
    void write_test(int ts, int n, Species& species);

    private:
    std::filesystem::path outputfolder;
    Domain& domain;
    std::ofstream file_data;
    std::ofstream file_ke;
    std::ofstream file_test;
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
