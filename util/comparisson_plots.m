filenames_o2 = {
    'mesh/solution_ShelfLife_O_2.txt', 
    'mesh/solution_Orchard_O_2.txt',
    'mesh/solution_Precooling_O_2.txt',
    'mesh/solution_Refrigerator_O_2.txt',
    'mesh/solution_OptimalCA_O_2.txt',
    'mesh/solution_DisorderInducing_O_2.txt'
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
    'mesh/solution_ShelfLife_CO_2.txt',
    'mesh/solution_Orchard_CO_2.txt'
    'mesh/solution_Precooling_CO_2.txt'
    'mesh/solution_Refrigerator_CO_2.txt'
    'mesh/solution_OptimalCA_CO_2.txt'
    'mesh/solution_DisorderInducing_CO_2.txt'
    };

%% Plot the solution 
for i = 1:length(filenames_co2)
    sol_o2 = readmatrix(filenames_o2{i}); 
    sol_co2 = readmatrix(filenames_co2{i});


    % graphic representation
    elements3   = Elements( : , 2:end ) ;
    coordinates = Nodes(:, 2:end); 
    figure('DefaultAxesFontSize',22)
    subplot(1, 2, 1)
    show(elements3,[],coordinates,full(sol_o2(1:end, 4)));
    title('oxygen')
    subplot(1, 2, 2)
    show(elements3,[],coordinates,full(sol_co2(1:end, 4)));
    title('carbon dioxide')

    saveas(gcf,strcat(states{i},'.png'))
end