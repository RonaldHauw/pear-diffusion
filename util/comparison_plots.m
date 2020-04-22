filenames_o2 = {
    'data/solutions/solution_ShelfLife_O_2.txt', 
    'data/solutions/solution_Orchard_O_2.txt',
    'data/solutions/solution_Precooling_O_2.txt',
    'data/solutions/solution_Refrigerator_O_2.txt',
    'data/solutions/solution_OptimalCA_O_2.txt',
    'data/solutions/solution_DisorderInducing_O_2.txt'
    };
states = {
    'ShelfLife', 
    'Orchard', 
    'Precooling', 
    'Refrigerator', 
    'OptimalCA', 
    'DisorderInducing'
    };
filenames_co2 = {
    'data/solutions/solution_ShelfLife_CO_2.txt',
    'data/solutions/solution_Orchard_CO_2.txt'
    'data/solutions/solution_Precooling_CO_2.txt'
    'data/solutions/solution_Refrigerator_CO_2.txt'
    'data/solutions/solution_OptimalCA_CO_2.txt'
    'data/solutions/solution_DisorderInducing_CO_2.txt'
    };



%% Plot the solution 
for i = 1:length(filenames_co2)
    
    run_software(states{i}); 
    %sol_o2 = readmatrix(filenames_o2{i}); 
    %sol_co2 = readmatrix(filenames_co2{i});
    %elements3   = Elements( : , 2:end ) ;
    %coordinates = Nodes(:, 2:end); 
    %[T_cel, n_u, n_v, name] = read_input(states{i});
    %show( [sol_o2(:, 4); sol_co2(:, 4)], name, elements3, coordinates, T_cel, n_u, n_v)

    % graphic representation
    
    %figure('DefaultAxesFontSize',22)
    %subplot(1, 2, 1)
    %show(elements3,[],coordinates,full(sol_o2(1:end, 4)));
    %title('oxygen')
    %subplot(1, 2, 2)
    %show(elements3,[],coordinates,full(sol_co2(1:end, 4)));
    %title('carbon dioxide')
    saveas(gcf,strcat('data/figs/',states{i},'.png'))
end