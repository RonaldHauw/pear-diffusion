
%% Select data points 
coarsest = 271
finest = 1571         % 
collect_data = 0    % should we collect data
residue_plot = finest-coarsest-1    % which level is used for plotting the residue 
coarse_grid = 10   % to which grid is everything interpolated (0 for a chosen set of points)


grid_finesses = [linspace(3, 140, 140-3+1), 150, 160, 170, 180, 190, 201, 251, 271, 379, 488];
condition = 'Precooling'; 

%% Make the grid and collect the data

addpath('util/')

if collect_data == 1

        %%% create the grid 
    for i = 1:length(grid_finesses)
        %pear_grid(grid_finesses(i)); 
        create_mesh( 'cc', 'pear', 0.025, grid_finesses(i))
        run_software(condition);
        copyfile(strcat('data/solutions/solution_',condition,'_O_2.txt'), strcat('data/solutions/solution_',condition,'_O_2_', int2str(grid_finesses(i)), '.txt')) 
        copyfile(strcat('data/solutions/solution_',condition,'_CO_2.txt'), strcat('data/solutions/solution_',condition,'_CO_2_', int2str(grid_finesses(i)), '.txt')) 

    end

end


%% Load solutions


if coarse_grid == 0
    coarse_sol_o2 = [0 0.00 0.03; 0 0.00 0.03]; 
    coarse_sol_co2 = [0 0.00 0.03; 0 0.00 0.03]; 
else 
    coarse_sol_o2  =readmatrix(strcat('data/solutions/solution_',condition,'_O_2_', int2str(coarse_grid),'.txt')); 
    coarse_sol_co2 = readmatrix(strcat('data/solutions/solution_',condition,'_CO_2_', int2str(coarse_grid),'.txt'));

end


sols_o2 = {}; 
sols_co2 = {}; 
for i = 1:length(grid_finesses)
    sol_o2 = readmatrix(strcat('data/solutions/solution_',condition,'_O_2_', int2str(grid_finesses(i)),'.txt')); 
    sol_co2 = readmatrix(strcat('data/solutions/solution_',condition,'_CO_2_', int2str(grid_finesses(i)),'.txt')); 
    sols_o2{i} = sol_o2; 
    sols_co2{i} = sol_co2; 
end


%% Interpolate to coarse grid

fine2coarse_sols_o2 = zeros(size(coarse_sol_o2, 1), length(sols_o2));
fine2coarse_sols_co2 = zeros(size(coarse_sol_co2, 1), length(sols_co2));

for i = 1:length(sols_o2)
    
    fine_sol_o2 = sols_o2{i}; 
    fine_sol_co2 = sols_co2{i}; 

    fine_interp_o2 = scatteredInterpolant(fine_sol_o2(:, 2),fine_sol_o2(:, 3),fine_sol_o2(:, 4));
    fine_interp_co2 = scatteredInterpolant(fine_sol_co2(:, 2),fine_sol_co2(:, 3),fine_sol_co2(:, 4));

    fine2coarse_sols_o2(:, i) = fine_interp_o2(coarse_sol_o2(:, 2), coarse_sol_o2(:, 3)); 
    fine2coarse_sols_co2(:,i) = fine_interp_co2(coarse_sol_co2(:, 2), coarse_sol_co2(:, 3)); 
end

%% calculate errors 

fine2coarse_errs_o2 = fine2coarse_sols_o2 - fine2coarse_sols_o2(:,end); 
fine2coarse_errs_co2 = fine2coarse_sols_co2 - fine2coarse_sols_co2(:,end);

o2_norms = vecnorm(fine2coarse_errs_o2)./ vecnorm(fine2coarse_sols_o2); % be careful for division with small numbers her!!
co2_norms = vecnorm(fine2coarse_errs_co2)./ vecnorm(fine2coarse_sols_co2);


%% 
radius = 0.01; 
max_elem_size = grid_finesses .^(-1); 
set(gca,'FontSize',24)
set(gca,'xscale','log')
set(gca,'yscale','log')
figure('DefaultAxesFontSize',22)
loglog(max_elem_size, o2_norms)
hold on 
scatter(max_elem_size, o2_norms)
title('Convergence of O2')
xlabel('Maximum element size')
ylabel('L2 norm of error')
figure('DefaultAxesFontSize',22)
loglog(max_elem_size, co2_norms)
hold on 
scatter(max_elem_size, co2_norms)
title('Convergence of CO2')
xlabel('Maximum element size')
ylabel('L2 norm of error')


figure, 
addpath('data/meshes')
create_mesh( 'cc', 'pear', 0.025, coarse_grid)
%pear_grid(coarse_grid)
load pear.mat
coordinates = Nodes(:, 2:3) ;
elements3   = Elements( : , 2:end ) ;
show_residual( [fine2coarse_errs_o2(:, residue_plot); fine2coarse_errs_co2(:, residue_plot)], condition, elements3, coordinates, 0, 0, 0)


