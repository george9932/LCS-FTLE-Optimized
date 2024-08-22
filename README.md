# LCS-FTLE-Optimized

LCS-FTLE-Optimized is a library for performing fast Lagrangian Coherent Structure (LCS) analysis of flow fields. It is a fork of the [lcs library](https://github.com/stevenliuyi/lcs) developed by Yi Liu. [OpenMP](http://www.openmp.org/) is supported for parallelization. 

It implements particle advection calculations for discrete flow velocity data, both forward and backward in time, and then uses the Finite-Time Lyapunov Exponents (FTLE) to obtain LCSs.  

A faster method that the classic FTLE calculations is used in this library. The unidirectional method described in the paper "[Fast computation of finite-time Lyapunov exponent fields for unsteady flows](https://doi.org/10.1063/1.3270044)" by Steven Brunton and Clarence Rowley is implemented here. The method ensures that unnecessary particle integrations are not repeated, leading to significant computational time savings.

## Development Environment

This project is set up to be run and debugged in Visual Studio Code (VSCode). It includes predefined tasks and launch configurations to facilitate convenient execution and debugging. Additionally, it is possible to run the project using the bash commands detailed later in this file.

### Requirements

- **Visual Studio Code**: Make sure you have [Visual Studio Code](https://code.visualstudio.com/) installed.
- **C/C++ Extension**: Install the [C/C++ extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) from Microsoft for code editing and debugging support.

### Tasks and Launch Configurations

The project includes tasks and launch configurations defined in the `.vscode` directory:

- **Tasks**: Used for building the project. You can run the tasks defined in the `tasks.json` file to compile the code.
- **Launch Configurations**: Used for debugging the application. The `launch.json` file contains configurations for running and debugging the application in VSCode.

Refer to the respective JSON files in the `.vscode` directory for further customization and configuration.

### Dependencies

Before compiling and running the program, ensure the following dependencies are installed:

1. **C++ Compiler**: The program requires a C++17 compatible compiler (e.g., `g++` version 5.0 or later) to compile the code.
   
2. **OpenMP**: The code uses [OpenMP](http://www.openmp.org/) for parallelization. Ensure your compiler supports OpenMP.
   
3. **Python 3**: Required for plotting and creating animations from the results.
   
4. **Python Libraries**: The following Python libraries are needed:
   - `matplotlib` for plotting FTLE results.
   - `numpy` for numerical operations.
   - `json`: Used for parsing and manipulating JSON data, particularly for reading and writing configuration files.
    - `os`: Provides a way to interact with the operating system, such as file and directory operations.
    - `timeit`: A utility for measuring the execution time of small code snippets to assess performance.
    - `multiprocessing`: Allows for parallel execution of tasks, enabling efficient utilization of multiple CPU cores.
    - `PIL`: A powerful image processing library used for opening, editing, and saving various image file formats.
    - `glob`: Facilitates file pattern matching, such as retrieving all files in a directory that match a certain pattern.
    - `sys`: Provides access to system-specific parameters and functions, often used for command-line arguments and exiting programs.

    The ones that are not part of the Python standard library can be installed via pip:
    
    ~~~~~~~~~~
    pip install matplotlib numpy Pillow
    ~~~~~~~~~~

5. **Git**: Required for cloning the repository.

## License

This project is licensed under the [GNU General Public License (GPL)](https://www.gnu.org/licenses/gpl-3.0.html). 

### Summary of the GNU GPL:

- You are free to use, modify, and distribute this software.
- Any modified versions must also be licensed under the GPL and made available to the public.
- You must provide a copy of the GPL license along with any distribution of the software.
- There is no warranty; the software is provided "as is" without any guarantees.

For detailed information, please refer to the `LICENSE` file in the repository.

## Installation
First, clone the repository from GitHub and access the repository folder:

~~~~~~~~~~
git clone https://github.com/george9932/LCS-FTLE-Optimized.git
cd LCS-FTLE-Optimized
~~~~~~~~~~

## Data Creation
In order to test the present library, it is possible to run some benchmark tests using analytical flow fields. If discrete data are available, this section can be skipped. The header file `src/velocity_function.hpp` includes 3 analytical velocity fields: the Bower model for meandering jet, the double-gyre model and the clockwise co-rotating double-gyre model. Discrete velocity data of these fields can be created by compiling and running `fast_computation/create_discrete_data.cpp` as:

~~~~~~~~~~
mkdir -p fast_computation/data
mkdir -p bin
g++ -g -Wall -O2 -std=c++17 -fopenmp fast_computation/create_discrete_data.cpp -o bin/CreateDiscreteData
bin/CreateDiscreteData
~~~~~~~~~~

The created data will be available in the `fast_computation/data` folder.

## JSON file
The discrete data are created and the calculations are done based on the values provided in the `fast_computation/sim_params.json` file that has the following structure (the present values are an example):

~~~~~~~~~~
{
    "x_min": 0.0,
    "x_max": 2.0,
    "y_min": 0.0,
    "y_max": 1.0,
    "nx": 500,
    "ny": 250,
    "data_nx": 500,
    "data_ny": 250,
    "t_min": 0,
    "t_max": 20,
    "data_delta_t": 0.2,
    "steps": 100,
    "file_prefix": "double_gyre_",
    "direction": "forward"
}
~~~~~~~~~~

The 2D rectangular domain is defined by `x_min`, `x_max`, `y_min`, `y_max`. The spatial discretization of the grid where the FTLE results are calculated is set by `nx` and `ny` . The spatial discretization of the velocity data is set with `data_nx` and `data_ny`. 

Velocity data from `t_min` to `t_max` with a time discretization of `data_delta_t` are created and written into files. Afterwards, the simulations run in t=[`t_min`,`t_max`]. The time direction for particle integration is specified by `direction` and can be either "forward" or "backward". The number of calculation steps is set by `steps`. 

All filenames start with the prefix specified in `file_prefix`. The data files can be found in the `data` folder. 

## Data File Structure
The data file structure follows a specific pattern. The first three lines of the file specify the number of elements in the x-direction `nx` (line 1), the number of elements in the y-direction `ny` (line 2), and the time associated with the saved field (line 3). Following these lines, the field data is stored sequentially.

For a 2D velocity field, the data corresponding to each grid point is saved as follows: the x-velocity is recorded on one line, followed by the y-velocity on the next line. The grid point at (`xmin`,`ymin`) is saved first.

Starting from this point, the y-coordinate is varied from `y_min` to `y_max` while keeping the x-coordinate fixed. Once all values for a given x-coordinate are recorded, the x-coordinate is incremented by dx, and the process repeats: the y-coordinate is varied again from `y_min` to `y_max`. This pattern continues until the data for the entire grid, up to (`x_max`,`y_max`), is saved.   

## Compilation and Run
Please make sure that C++17 standard is supported by the compiler, and then compile the program:

~~~~~~~~~~
make
~~~~~~~~~~

The compiled program will be found in the `bin` folder. The program can be run as:

~~~~~~~~~~
bin/DiscreteFastComputation
~~~~~~~~~~

The primary output of the program is the FTLE results, located in `fast_computation/results/ftle`. Additionally, the program calculates and stores step flow maps as binary files in the `fast_computation/step_flow_maps` directory. These maps represent the positions of fluid particles after a single time step, originating from a uniform grid.

## Plot and Animation
The results can be plotted by:

~~~~~~~~~~
python3 fast_computation/results/plot_ftle.py
~~~~~~~~~~

The contour levels and the colormap can be defined inside the main function of `fast_computation/results/plot_ftle.py`. The plots will then be available in the `fast_computation/results/ftle` folder.

Afterwards, an animation can be created from the plots by:

~~~~~~~~~~
python3 fast_computation/results/animation_ftle.py
~~~~~~~~~~

The number of Frames Per Second (FPS) can be specified in the main function of `fast_computation/results/animation_ftle.py`. The animation will then be available in the `fast_computation/results/ftle` folder.

Make sure that the values defined in `fast_computation/sim_params.json` are the same as in the simulations.
