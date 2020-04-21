# The Pear Project (2019-2020)
## Solving a macroscale respiration–diffusion model

### Authors
- Agathe, van Lamsweerde;
- Octave, Oliviers;
- Ronald, Hauwaerts.

Our code is free to use and modify. 

### Installation
#### Needed software:
- A C++ compiler which works with the C++17 standard such as *g++* (https://gcc.gnu.org);
- Matlab (https://www.mathworks.com/products/matlab.html).
#### Needed libraries:
- The C++ *standard library* which is usually preïnstalled;
- A copy of *eigen* (http://eigen.tuxfamily.org) at the location "software/". 
    
    
        pear_group11$ cd software 
        software$ git clone https://gitlab.com/libeigen/eigen.git 
          
#### Compilation: 
##### 1. Shell script 
Open a command line window in the root directory and use the following command. 

    pear_group11$ ./util/compile.sh 
The script will use the flags -03 and -Wall. The compilation can use a lot of memory, if you have Cmake installed, the following method may be preferred.  
##### 2. Via Cmake  
Open a command line window in the root directory and use the following commands. 

    pear_group11$ cmake . 
    pear_group11$ make 
##### 3. Manually 
Open a command line window in the root directory and use the following command. 

    pear_group11$ g++ -Wall -std=c++17 -O3    software/component.hpp \
                                   software/diffusion.hpp \
                                   software/grid.hpp \
                                   software/rdc.hpp \
                                   software/reaction.hpp \
                                   software/nlsolver.hpp \
                                   software/main.cpp
        
### Running the software
#### A first example
The Matlab script "run_software.m" runs the executable and automatically provides some plots of the results. 
Open the Matlab command line window in in the root directory: 

    >> run_software("Orchard")
    
This will run a simulation for the "Orchard" configuration such as described in the assignment. The other configurations can be 
tested by replacing "Orchard" with the appropriate name. 

#### Creating a custom mesh 
The matlab function "creat_mesh.m" creates a mesh and stores it in the right location for it to be read by the software executable.
The mesh which is used by the C++ code is stored in data/meshes and always has the name 'pear'. 

    >> addpath("util")
    >> create_mesh('pear', 'pear', 0.025, 30 ) 
    
The last line creates a 'pear'-shaped grid which is named 'pear', has a radius of 0.025 meters and uses a node distance of 
approximately 0.025/30. Other predefined shapes are 

- double circle; 
- circle triangle;
- circle.

#### Manually calling the executable
To run the executable with it's default hyperparameters use:

    pear_group11$ ./pear_diffusion 

Some useful hyperparameters are for example: 

    pear_group11$ ./pear_diffusion -ShelfLife -res_pred 5e-15 -res_new 1e-17
   
This sets the shelf life parameters and sets the residual tubes for the non-linear solver, effictively bounding the residual between 
5e-15 in the prediction step and stopping the Newton iterations if the residual is smaller than 1e-17. 

### Expanding the code 

Code is constantly in flux and we've made our best to make the source code easy to interpret and change. For example 
    
   - the code has a logical and clear object oriented structure;
   - we make use of abstraction as much as possible to make the code reusable;   
   - the code is templated such that the data type can be changed very quick. 

#### Example: adding a component 
Let's say we want to model not only O2 and CO2 but also ethylene. This example provides some guidance on how this can be done
In main.cpp set the number of components

    int nb_comp = 3; 
    
and declare a new 'component'

    pear::component<d_type, vec_type, mat_type> ethylene("Ethylene", grid, conc, 1, 2);

The last two numbers indicate that we use no stride (stride=1) and that ethylene is the third (2+1) component in the solution vector. 
Fill in the parameters with which etylene will diffuse and declare the 'diffusion' equations:

    std::vector<d_type> diffusion_ethylene_param = {sigma_e_r, sigma_e_z, r_v, c_e_amb};
    pear::diffusion<d_type, vec_type, mat_type> diff_ethylene(ethylene, grid, diffusion_co2_param);

Write a (small) class similar to 'respiration' which provides the reaction kinetics and their partial derivatives. 
Declare the 'reaction' equations: 

    pear::reaction<d_type , vec_type, mat_type> react_o2_ethylene(co2, o2, grid, kinetics_o2_ethylene);

Currently the 'reaction' class only supports a reaction between components but this could be easily extended. Modify the 
'reaction-diffucion-convection' class which makes abstraction of the underlying diffusion and reaction and behaves as a function. 

    pear::rdc<d_type, vec_type, mat_type> equation(diff_o2, diff_co2, diff_ethylene, react_co2_o2, react_o2_ethylene);

This short example shows that none of the classes change significantly. For example, we could add a new component which only diffused in 10 minutes. The only two 'large' changes are:

1. Writing a class which contains the specific reaction kinetics;
2. Modifying the 'reaction-diffusion-convection' to include the new diffusion and reaction. 



