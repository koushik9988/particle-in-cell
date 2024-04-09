#include "output.h"
#include <memory>


Output::Output(const std::filesystem::path& outputfolder, Domain& domain) : outputfolder(outputfolder), domain(domain) 
{

    std::filesystem::remove_all(outputfolder);  // Clear previous output
    std::filesystem::create_directories(outputfolder);
    //std::filesystem::create_directories(outputfolder / "files");

    //file_data.open(outputfolder / "files" / "results_" + std::to_string(domain.alpha) + ".txt");
    file = H5File(outputfolder / "result.h5", H5F_ACC_TRUNC);
    //file = H5File(outputfolder / "result.h5",H5F_ACC_RDWR);
    
    
    if (file.getId() < 0) 
    {
        throw std::runtime_error("Error opening HDF5 file");
    }
    //create groups
    //particle_group1 = file.createGroup("/electron");
    //particle_group2 = file.createGroup("/ion");
    field_data_group = file.createGroup("/fielddata");
    time_group = file.createGroup("/time_var");
    metadata_group = file.createGroup("/metadata");


    
    //-----------
    int sp_no = domain.species_no;
    int t_step = int(domain.NUM_TS/domain.write_interval) + 1;
    
    store_ke = new double*[t_step];
    for (int i=0; i < t_step; i++)
    {
        store_ke[i] = new double[sp_no];
    }
    //clear memory of all pointers
    for (int i = 0; i < t_step; i++)
    {
        for (int j = 0; j < sp_no ; j++)
        {
            store_ke[i][j]= 0;  
        }
    }
    //--------------------------------
}

//-------
void Output::write_metadata(int NC, int NUM_TS, int write_int, int write_int_phase, double DT, int nE, int nI, int nN, int nB, double Te, double Tm, double Tb, double alpha, double beta, double mI, double mN, double mB, double density)
{
    // Write metadata attributes within the metadata group
    //Group metadata_group = file.openGroup("/metadata");

    metadata_group.createAttribute("NC", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &NC);
    metadata_group.createAttribute("NUM_TS", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &NUM_TS);
    metadata_group.createAttribute("write_int", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &write_int);
    metadata_group.createAttribute("write_int_phase", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &write_int_phase);
    metadata_group.createAttribute("DT_coeff", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &DT);
    metadata_group.createAttribute("nE", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &nE);
    metadata_group.createAttribute("nI", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &nI);
    metadata_group.createAttribute("nN", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &nN);
    metadata_group.createAttribute("nB", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &nB);
    metadata_group.createAttribute("Te", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &Te);
    metadata_group.createAttribute("Tm", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &Tm);
    metadata_group.createAttribute("Tb", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &Tb);
    metadata_group.createAttribute("alpha", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &alpha);
    metadata_group.createAttribute("beta", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &beta);
    metadata_group.createAttribute("mI", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &mI);
    metadata_group.createAttribute("mN", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &mN);
    metadata_group.createAttribute("mB", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &mB);
    metadata_group.createAttribute("density", PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_DOUBLE, &density);

    metadata_group.close();
}

//---------

void Output::write_field_data(int ts)
{
    
    std::string pot_group_name  = "pot";
    std::string efield_group_name  = "efield";

    Group pot_subgroup;
    Group efield_subgroup;

    if (!field_data_group.exists(pot_group_name))
    {
        // Subgroup doesn't exist, create it
        pot_subgroup = field_data_group.createGroup(pot_group_name);
        efield_subgroup = field_data_group.createGroup(efield_group_name);
    }
    else
    {
        // Subgroup already exists, retrieve it
        pot_subgroup = field_data_group.openGroup(pot_group_name);
        efield_subgroup = field_data_group.openGroup(efield_group_name);
    }
    
    //pot_subgroup = field_data_group.createGroup(subgroup_name);
    

    hsize_t ni = domain.ni;

    hsize_t dims_den[1] = {ni};

    hsize_t Rank = 1;

    DataSpace dataspace_pot(Rank, dims_den);
    DataSpace dataspace_efield(Rank, dims_den);

    H5::DataSet dataset_pot = pot_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_pot);
    H5::DataSet dataset_efield = efield_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_efield);

    // Prepare data buffer
    std::vector<double> pot(ni);
    std::vector<double> efield(ni);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int i = 0; i < ni; ++i) 
    {
        pot[i] = domain.phi[i];
        efield[i] = domain.ef[i];  
    }

    // Write the density data to the dataset
    dataset_pot.write(pot.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_efield.write(efield.data(), H5::PredType::NATIVE_DOUBLE);
}

//-----------------den data----------------------------------------
void Output::write_den_data(int ts,  Species& species)
{
    
    std::string subgroup_name = "den_" + species.name;

    Group den_subgroup;

    // Check if the group already exists in the map
    auto it = den_subgroups.find(species.name);
    if (it != den_subgroups.end()) 
    {
        // Group already exists, retrieve it from the map
        den_subgroup = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        den_subgroup = field_data_group.createGroup(subgroup_name);
        den_subgroups[species.name] = den_subgroup;
    }
    
    hsize_t ni = domain.ni;

    hsize_t dims_den[1] = {ni};

    hsize_t Rank = 1;

    DataSpace dataspace_den(Rank, dims_den);

    H5::DataSet dataset_den = den_subgroup.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_den);

    // Prepare data buffer
    std::vector<double> density(ni);

    //density = species.den;

    // Flatten 2D array into a 1D vector for writing into the dataset
    for (int i = 0; i < ni; ++i) 
    {
        density[i] = species.den[i];
    }

    // Write the density data to the dataset
    dataset_den.write(density.data(), H5::PredType::NATIVE_DOUBLE);
}

void Output::write_particle_data(int ts, Species& species)
{
    std::string group_name = "particle_" + species.name;
    Group particle_group;

    // Check if the group already exists in the map
    auto it = particle_groups.find(species.name);
    if (it != particle_groups.end()) 
    {
        // Group already exists, retrieve it from the map
        particle_group = it->second;
    }
    else 
    {
        // Group doesn't exist, create it and store it in the map
        particle_group = file.createGroup(group_name);
        particle_groups[species.name] = particle_group;
    }

    // Create datasets for particle positions and velocities
    hsize_t dims_phase[2] = {species.part_list.size(), 2}; //here 2 => 1 for pos and 1 for vel
    //hsize_t dims_vel[2] = {species.part_list.size(), 1};
    hsize_t Rank = 2;

    DataSpace dataspace_phase(Rank, dims_phase);
    //DataSpace dataspace_vel(Rank, dims_vel);

    H5::DataSet dataset_phase = particle_group.createDataSet(std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_phase);
    //H5::DataSet dataset_vel = particle_group.createDataSet("vel" + std::to_string(ts), H5::PredType::NATIVE_DOUBLE, dataspace_vel);

    // Write particle positions and velocities to the datasets
    std::vector<double> phase_data;
    //std::vector<double> velocities;

    for (Particle& p : species.part_list) 
    {
        phase_data.push_back(p.pos);
        phase_data.push_back(p.vel);
    }

    dataset_phase.write(phase_data.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vel.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
}

void Output::write_ke()
{
    // Define the name for the dataset
    std::string datasetName = "kinetic_energy";

    // Define the dimensions of the dataset
    hsize_t ny = int(domain.NUM_TS/domain.write_interval) + 1; //rows
    hsize_t nx = 3; //column

    hsize_t dims_energy[2] = {ny, nx}; // Assuming all rows have the same length

    // Create dataspace for the dataset
    DataSpace dataspace_energy(2, dims_energy);

    // Create the dataset within the time_group
    H5::DataSet dataset_energy = time_group.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace_energy);

    // Allocate memory for the data (outside the loop for efficiency)
    double* ke_data = new double[nx * ny];

    // Fill the data array
    for (hsize_t i = 0; i < ny; ++i)
    {
        for (hsize_t j = 0; j < nx; ++j)
        {
            ke_data[i * nx + j] = store_ke[i][j];
        }
    } 

    // Write the data to the dataset
    dataset_energy.write(ke_data, H5::PredType::NATIVE_DOUBLE);

    // Deallocate memory after writing
    //delete[] ke_data;
}
//*/
void Output::printmatrix(int row, int col, double **matrix)
{
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

