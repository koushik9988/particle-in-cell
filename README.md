# Electrostatic 1D Particle-in-Cell (PIC) Code

This repository contains an electrostatic 1D Particle-in-Cell (PIC) code developed for simulating plasma systems. The code is capable of simulating basic plasma phenomena.

## Requirements
- Python3
- GNU C++ compiler(g++)
- GNU make
- git
- Matplotlib
- NumPy
- Scipy


### Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/koushik9988/particle-in-cell.git
    ```

2. Navigate to the directory:
    ```bash
    cd your_repository
    ```

3. Build the code using make:
    ```bash
    make clean
    ```
    ```bash
    make all
    ```

### Running the Code
1. Configure the simulation parameters in the `input.ini` file.
2. Run the code:
    ```bash
    g++ ./pic
    ```

# Explanation of input.ini file Parameters

The input file `input.ini` contains parameters for configuring the simulation. Each section corresponds to a different aspect of the simulation setup.

## `[file]`

- **output**: Specifies the directory where the output data will be stored.

## `[time]`

- **NUM_TS**: Total number of time steps for the simulation.

## `[diagnostics]`

- **write_interval**: Interval for writing density and field data in result.txt file.
- **write_interval_phase**: Interval for writing phase-space data file.
- **write_diagnostics**: Interval for writing diagnostic outputs.
- **DT_coeff**: Coefficient for the time step.
- **write_flag**: Flag for controlling data writing: 
  - 1: Write both phase and field data.
  - 2: Write only field data.
  - 3: Write only phase data.
  - 0: Write no data.

## `[domain]`

- **NC**: Number of cells in the domain.
- **x0**: Initial position of the domain.

## `[population]`

- **nParticlesE**: Number of electrons loaded into the domain.
- **nParticlesI**: Number of ions loaded into the domain.
- **nParticlesN**: Number of negative ions loaded into the domain.
- **nParticlesB**: Number of background particles.
- **tempE**: Temperature of electrons.
- **tempI**: Temperature of ions.
- **tempN**: Temperature of negative ions.
- **tempB**: Temperature of background particles.
- **massE**: Mass of electrons.
- **massI**: Mass of ions.
- **massN**: Mass of negative ions.
- **massB**: Mass of background particles.

## `[simulation]`

- **v_i**: Ion streaming velocity.
- **v_e**: Electron streaming velocity.
- **v_b**: Beam particle streaming velocity.
- **v_n**: Negative ion streaming velocity.
- **density**: Plasma density.
- **alpha**: Fraction of negative ion to background positive ion
- **beta**: fraction of negative beam to background negative ion 
- **bc**: Boundary condition.
   - 1: pbc for periodic boundary.
   - 2: open for open boundary condition.



![dispersion](https://github.com/koushik9988/particle-in-cell/assets/55924787/5d278d78-2755-4293-bf18-4f8a09789b8c)
