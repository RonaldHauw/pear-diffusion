%% Group 11 - April 13th 2020
% Finite element method for two-dimensional reaction-diffusion system.

%% Compile the C++ code if executable does not exist
if ~ isfile('./pear_diffusion')
    ! ./compile.sh
end

%% Initialise the mesh grid
run('matlab/hct_grid.m')

%% Solve using C++
! ./pear_diffusion -maxit 100 -ShelfLife


%% Plot the solution 
sol_o2 = readmatrix('data/solutions/solution_ShelfLife_O_2.txt'); 
sol_co2 = readmatrix('data/solutions/solution_ShelfLife_CO_2.txt');

addpath('matlab/') 

% graphic representation
elements3   = Elements( : , 2:end ) ;
coordinates = Nodes(:, 2:end); 

show( [sol_o2(:, 4); sol_co2(:, 4)], 'shelf life', elements3, coordinates, [0, 0] )