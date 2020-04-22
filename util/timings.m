%%% create the grid 
grid_finesses = [10, 20, 30, 50, 100, 150, 200, 400, 600, 100];
condition = 'Refrigerator';
nb_samples = 2; 

c_times = zeros(length(grid_finesses), 1);
m_times = zeros(length(grid_finesses), 1);
addpath("prototype")

for i = 1:length(grid_finesses)
    create_mesh( 'pear', 'pear', 0.025, grid_finesses(i)) % create grid 
    run_software(condition);
    c_t = 0.0; 
    for j = 1:nb_samples
        tic(); 
        run_software(condition);
        c_t = c_t+ toc(); 
    end
    c_times(i) = c_t/nb_samples; 
    m_t = 0.0; 
    for j = 1:nb_samples
        tic(); 
        run(condition); 
        m_t = m_t + toc(); 
    end
    m_times(i) = m_t/nb_samples; 
    
    copyfile(strcat('data/solutions/solution_',condition,'_O_2.txt'), strcat('data/solutions/solution_',condition,'_O_2_', int2str(grid_finesses(i)), '.txt')) 
    copyfile(strcat('data/solutions/solution_',condition,'_CO_2.txt'), strcat('data/solutions/solution_',condition,'_CO_2_', int2str(grid_finesses(i)), '.txt')) 

end