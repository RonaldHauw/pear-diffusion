%% Group 11 - April 13th 2020
%
% Solve reaction-diffusion system with finite element method in C++.
%
% Either run simulation by name as
%     >> run( name )
% where 'name' is the string 'orchard', 'shelf life', 'refrigerator', 
% 'precooling', 'disorder inducing' or 'optimal CA',
%
% or run default 'refrigerator' simulation ( T_cel = 7, n_u = 0.208, n_v = 0 )
%    >> run

function run_software( varargin )
    clc

    %% Read input
    addpath('util/')
    [T_cel, n_u, n_v, name] = read_input( varargin{:} ) ;

    
    %% Compile the C++ code if executable does not exist
    if ~ isfile('pear_diffusion')
        ! ./util/compile.sh
    end
    

    %% Solve using C++
    if contains(lower(name), 'orchard')
        sim = 'Orchard' ;
    elseif contains(lower(name), 'shelf')
        sim = 'ShelfLife' ;
    elseif contains(lower(name), 'refr')
        sim = 'Refrigerator' ;
    elseif contains(lower(name), 'pre')
        sim = 'Precooling' ;
    elseif contains(lower(name), 'diso')
        sim = 'DisorderInducing' ;
    elseif contains(lower(name), 'optim')
        sim = 'OptimalCA' ;
    end
    
    command = strcat('./pear_diffusion', ' -maxit 100', ' -', sim, ' -res_pred 5e-15', ' -res_new 1e-16');
    system(command);

    %% Plot the solution
    path = 'data/solutions/solution_' ;
    
    sol_o2   = readmatrix( strcat( path, sim, "_O_2.txt" ));
    sol_co2  = readmatrix( strcat( path, sim, "_CO_2.txt" ));

    % graphic representation
    addpath('data/meshes')
    load pear.mat
    coordinates = Nodes(:, 2:3) ;
    elements3   = Elements( : , 2:end ) ;

    show( [sol_o2(:, 4); sol_co2(:, 4)], name, elements3, coordinates, T_cel, n_u, n_v)
end