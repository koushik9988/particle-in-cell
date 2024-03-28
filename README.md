# Electrostatic 1D Particle-in-Cell (PIC) Code

This repository contains an electrostatic 1D Particle-in-Cell (PIC) code developed for simulating plasma systems. The code is capable of simulating basic plasma phenomena.

## Requirements
- Python3
- Matplotlib
- NumPy


### Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/your_username/your_repository.git
    ```

2. Navigate to the directory:
    ```bash
    cd your_repository
    ```

3. Build the code using make:
    ```bash
    make
    ```

### Running the Code
1. Configure the simulation parameters in the `input.ini` file.
2. Run the code:
    ```bash
    python main.py
    ```

### Input Configuration (input.ini)
The `input.ini` file contains parameters for configuring the simulation. Parameters include:
- **Simulation duration**
- **Number of particles**
- **Grid size**
- **Time step**
- **Initial conditions**

Refer to the comments in the `input.ini` file for detailed information on each parameter.

## Benchmarking
### Langmuir Wave Dispersion
To benchmark with the Langmuir dispersion relation, run the code with appropriate parameters and compare the simulation results with the theoretical dispersion relation.

### Two-Stream Instability
Visualize the stability phase plot to observe the behavior of the two-stream instability.

## Examples
Include example plots and results from simulations.

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License
This project is licensed under the [MIT License](LICENSE).

![dispersion](https://github.com/koushik9988/particle-in-cell/assets/55924787/5d278d78-2755-4293-bf18-4f8a09789b8c)
