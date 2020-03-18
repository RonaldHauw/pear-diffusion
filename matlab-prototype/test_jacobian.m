%% Group 11 - March 12th 2020
% test correctness of Jacobian with finite difference approximation

clear all

% load variables of system
load('workspace') ;
% choose simple circular domain
load meshes/HalfCircleMesh.mat
load meshes/HalfCircleMesh_Data.mat
% shrink size of mesh
coordinates = mesh.Nodes'/15 ;
elements3   = Elements( : , 2:end ) ;
% number of vertices
M           = size(coordinates, 1) ;

% create random vector of concentrations
C   = 100*randn( 2*M, 1 ) ;
% evaluate Jacobian in the random vector
J_c = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;

% maximal exponent
E = 16 ; path = zeros(1, E) ;
% store max difference between jacobian and approximation
for e = 1:E
    % perturbation
    %d_C   = randn(size(C))*10^(-e) ;
    d_C   = 10^(-e) ;

    H_mdc = assemble_H( coordinates, elements3, C - d_C, R_u, R_v ) ;
    H_pdc = assemble_H( coordinates, elements3, C + d_C, R_u, R_v ) ;

    path(e) = max(max( (H_pdc-H_mdc)./(2*d_C') - J_c )) ;
end

% plot convergence
figure
box on
hold on
loglog(10.^(-1:-1:-E), path, 'k-') ;
loglog(10.^(-1:-2:-E), path(end) + 10.^(-1:-2:-E), 'b--')
hold off
legend('Empirical convergence', 'Expected convergence')
title('Convergence of Finite difference approximaiton of the Jacobian')
xlabel('Order of perturbation on C')
ylabel({'Difference between', 'Jacobian and approximation'})