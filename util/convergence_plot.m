%% Load solutions

sol_o2_5 = readmatrix('mesh/solution_O_2_5.txt'); 
sol_co2_5 = readmatrix('mesh/solution_CO_2_5.txt');

sol_o2_10 = readmatrix('mesh/solution_O_2_10.txt'); 
sol_co2_10 = readmatrix('mesh/solution_CO_2_10.txt');

sol_o2_15 = readmatrix('mesh/solution_O_2_15.txt'); 
sol_co2_15 = readmatrix('mesh/solution_CO_2_15.txt');

sol_o2_20 = readmatrix('mesh/solution_O_2_20.txt'); 
sol_co2_20 = readmatrix('mesh/solution_CO_2_20.txt');

sol_o2_25 = readmatrix('mesh/solution_O_2_25.txt'); 
sol_co2_25 = readmatrix('mesh/solution_CO_2_25.txt');

sol_o2_30 = readmatrix('mesh/solution_O_2_30.txt'); 
sol_co2_30 = readmatrix('mesh/solution_CO_2_30.txt');

sols_o2 = {sol_o2_5, sol_o2_10, sol_o2_15, sol_o2_20, sol_o2_25, sol_o2_30}; 
sols_co2 = {sol_co2_5, sol_co2_10, sol_co2_15, sol_co2_20, sol_co2_25, sol_co2_30}; 

%% Convert all to coarsest

coarse_sol_o2 = sol_o2_5; 
coarse_sol_co2 = sol_co2_5; 

fine2coarse_sols_o2 = zeros(size(sol_o2_5, 1), length(sols_o2));
fine2coarse_sols_co2 = zeros(size(sol_co2_5, 1), length(sols_co2));

for i = 1:length(sols_o2)
    
    fine_sol_o2 = sols_o2{i}; 
    fine_sol_co2 = sols_co2{i}; 

    fine_interp_o2 = scatteredInterpolant(fine_sol_o2(:, 2),fine_sol_o2(:, 3),fine_sol_o2(:, 4));
    fine_interp_co2 = scatteredInterpolant(fine_sol_co2(:, 2),fine_sol_co2(:, 3),fine_sol_co2(:, 4));

    fine2coarse_sols_o2(:, i) = fine_interp_o2(sol_o2_5(:, 2), sol_o2_5(:, 3)); 
    fine2coarse_sols_co2(:,i) = fine_interp_co2(sol_co2_5(:, 2), sol_co2_5(:, 3)); 
end

%% calculate errors 

fine2coarse_errs_o2 = fine2coarse_sols_o2 - fine2coarse_sols_o2(:,end); 
fine2coarse_errs_co2 = fine2coarse_sols_co2 - fine2coarse_sols_co2(:,end);

o2_norms = vecnorm(fine2coarse_errs_o2);
co2_norms = vecnorm(fine2coarse_errs_co2);


%% 
radius = 0.02; 
grid_finesses = [1./5 1./10 1./15 1./20 1./25 1./30];
max_elem_size = grid_finesses*1.5*radius; 

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




