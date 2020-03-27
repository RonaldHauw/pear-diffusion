%% Compile the C++ code
%! ./compile.sh

%% Initialise the mesh grid
run('matlab/hct_grid.m')

%% Solve using C++

! ./pear_diffusion -maxit 100   -ShelfLife


%% Plot the solution 
sol_o2 = readmatrix('data/solutions/solution_ShelfLife_O_2.txt'); 
sol_co2 = readmatrix('data/solutions/solution_ShelfLife_CO_2.txt');

addpath('matlab/') 

% graphic representation
elements3   = Elements( : , 2:end ) ;
coordinates = Nodes(:, 2:end); 
figure()
subplot(1, 2, 1)
show(elements3,[],coordinates,full(sol_o2(1:end, 4)));
title('oxygen')
subplot(1, 2, 2)
show(elements3,[],coordinates,full(sol_co2(1:end, 4)));
title('carbon dioxide')
