#include "output.h"

Output::Output(const std::filesystem::path& outputfolder, Domain& domain) :
    outputfolder(outputfolder), domain(domain) {

    std::filesystem::remove_all(outputfolder);  // Clear previous output
    std::filesystem::create_directories(outputfolder);
    std::filesystem::create_directories(outputfolder / "files");

    //file_data.open(outputfolder / "files" / "results_" + std::to_string(domain.alpha) + ".txt");
    file_data.open(outputfolder / "files" / ("Results.txt"));

    if (!file_data.is_open()) {
        throw std::runtime_error("Error opening file_data");
    }

    //file_ke.open(outputfolder / "files" / "ke_" + std::to_string(domain.alpha) + ".txt");
    file_ke.open(outputfolder / "files" / ("ke_" + std::to_string(domain.alpha) + ".txt"));

    if (!file_ke.is_open()) {
        throw std::runtime_error("Error opening file_ke");
    }
}

void Output::write_particle(int ts, Species& species) 
{   
    std::filesystem::path outputFilepath;

    std::string name = species.name;
    
    outputFilepath = outputfolder / (name[0] + std::to_string(ts) + ".txt");
    
    std::ofstream outputFile(outputFilepath);

    std::ostringstream buffer;

    // Check if the file is opened successfully
    if (!outputFile.is_open()) 
    {
        std::cerr << "Error opening file: " << outputFilepath << std::endl;
        return; // Exit with an error
    }

    for (const Particle& p : species.part_list) 
    {
        //outputFile << p.pos << "\t" << p.vel << std::endl;
        buffer << p.pos << "\t" << p.vel << std::endl;
    }

    outputFile<<buffer.str();
}

void Output::write_data(int ts, std::vector<Species>& species_list) 
{
    // Print data for each grid point in a columnar format
    for (int i = 0; i < domain.ni; i++) 
    {
        file_data << int(i * domain.dx);  // Print the grid point value

        for (auto &sp : species_list) 
        {
            file_data << "\t\t" << std::fixed << std::setprecision(5) << sp.den[i];  // Print density with 2 decimal places
        }

        file_data << "\t\t" << std::fixed << std::setprecision(5) << domain.phi[i] << "\t\t" << std::fixed << std::setprecision(5) << domain.ef[i];  // Print electric field and potential with 2 decimal places
        file_data << std::endl;  // Move to the next line after each grid point
    }

    file_data << std::endl;  // Add an extra empty line at the end
}

void Output::write_ke(int ts,std::vector<Species> &species_list)
{
    std::vector<Species>::iterator it;

    it = std::find_if(species_list.begin(),species_list.end(), [](Species &sp) { return sp.name == "electron";});
    /*
    Syntax
    InputIterator find_if (InputIterator first, InputIterator last, UnaryPredicate pred);
    Parameters:-
    first, last: range which contains all the elements between first and last, including the element pointed by first but not the element pointed by last.
    pred: Unary function that accepts an element in the range as an argument and returns a value in boolean.
    Return Value
    This function returns an iterator to the first element in the range [first, last) for which pred(function) returns true. If no such element is found, the function returns last.
    */

    if (it != species_list.end())
    {
        Species &electron = *it;

        //for(Species &sp : species_list)
        //{
            //file_ke<<sp.name<<"\t \t";
        //}
         //file_ke<<std::endl;

        file_ke<<ts*domain.DT<<"\t"<<domain.ComputePE(electron);
        
        for(Species &sp : species_list)
        {
            file_ke<<"\t \t"<<std::fixed << std::setprecision(5)<<sp.Compute_KE(electron);
        }
        file_ke<<std::endl;
    }
}

