#include <vector>
#include "Species.h" // Assuming Species class definition is included here
#include "H5Cpp.h"   // Assuming HDF5 library is included here

class Output {
private:
    std::vector<double> time_steps_accumulated;
    std::vector<double> ke_data_accumulated;
    H5::Group time_group; // Assuming time_group is defined here
    // Add other necessary member variables
    
public:
    // Other member functions
    
    void accumulate_and_write_ke_data(int NUM_TS, int write_interval, const Domain& domain, std::vector<Species>& species_list) {
        // Accumulate kinetic energy data during the main loop
        for (int ts = 0; ts < NUM_TS + 1; ts++) {
            // Your main loop code here...
            
            if (ts % write_interval == 0) {
                accumulate_ke_data(ts, domain, species_list);
            }
        }
        
        // Write accumulated data to the dataset
        write_accumulated_ke_data(domain, species_list);
    }

private:
    void accumulate_ke_data(int ts, const Domain& domain, const std::vector<Species>& species_list) {
        // Find the electron species in the species list
        auto it = std::find_if(species_list.begin(), species_list.end(), [](const Species& sp) { return sp.name == "electron"; });

        if (it != species_list.end()) {
            // Electron species found
            const Species& electron = *it;

            // Prepare a temporary vector to store data for one timestep
            std::vector<double> ke_data_single_timestep(species_list.size() + 1); // Allocate space for time and KE data

            // Store time step in the first column
            ke_data_single_timestep[0] = ts * domain.DT;

            // Calculate and store kinetic energy for each species
            for (size_t i = 0; i < species_list.size(); ++i) {
                ke_data_single_timestep[i + 1] = species_list[i].Compute_KE(electron);
            }

            // Accumulate data for this timestep
            time_steps_accumulated.push_back(ke_data_single_timestep[0]);
            ke_data_accumulated.insert(ke_data_accumulated.end(), ke_data_single_timestep.begin() + 1, ke_data_single_timestep.end());
        }
    }

    void write_accumulated_ke_data(const Domain& domain, const std::vector<Species>& species_list) {
        if (time_steps_accumulated.empty() || ke_data_accumulated.empty()) {
            // Nothing to write
            return;
        }

        // Prepare dataset dimensions
        hsize_t dims_ke[2] = { static_cast<hsize_t>(time_steps_accumulated.size()), static_cast<hsize_t>(species_list.size() + 1) };
        hsize_t rank = 2;
        DataSpace dataspace_ke(rank, dims_ke);

        // Open the dataset (create it if it doesn't exist)
        H5::DataSet dataset_ke;
        if (time_group.exists("ke_energy")) {
            dataset_ke = time_group.openDataSet("ke_energy");
        } else {
            dataset_ke = time_group.createDataSet("ke_energy", H5::PredType::NATIVE_DOUBLE, dataspace_ke);
        }

        // Combine time steps and kinetic energy data
        std::vector<double> all_data;
        all_data.reserve(time_steps_accumulated.size() + ke_data_accumulated.size());
        all_data.insert(all_data.end(), time_steps_accumulated.begin(), time_steps_accumulated.end());
        all_data.insert(all_data.end(), ke_data_accumulated.begin(), ke_data_accumulated.end());

        // Write all data to the dataset
        dataset_ke.write(all_data.data(), H5::PredType::NATIVE_DOUBLE);

        // Clear accumulated data for the next simulation
        time_steps_accumulated.clear();
        ke_data_accumulated.clear();
    }
};

