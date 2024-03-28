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
    make all
    ```

### Running the Code
1. Configure the simulation parameters in the `input.ini` file.
2. Run the code:
    ```bash
    g++ ./pic
    ```

### Input Configuration (input.ini)
The `input.ini` file contains parameters for configuring the simulation. Parameters include:
- **Simulation duration**
- **Number of particles**
- **Grid size**
- **Time step**
- **Initial conditions**

Refer to the comments in the `input.ini` file for detailed information on each parameter.


![dispersion](https://github.com/koushik9988/particle-in-cell/assets/55924787/5d278d78-2755-4293-bf18-4f8a09789b8c)
