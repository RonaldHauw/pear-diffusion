% @Author: OctaveOliviers
% @Date:   2020-03-18 20:55:56
% @Last Modified by:   OctaveOliviers
% @Last Modified time: 2020-03-20 20:23:22

%% Group 11 - March 12th 2020
% test correctness of Jacobian with finite difference approximation

clear all

% load variables of system
load('workspace') ;
% choose simple circular domain
load meshes/HalfCircleMesh.mat
load meshes/HalfCircleMesh_Data.mat
% shrink size of mesh
coordinates = (mesh.Nodes'/15) ;
elements3   = Elements( : , 2:end ) ;
% number of vertices
M           = size(coordinates, 1) ;

% create random vector of concentrations
C   = randn( 2*M, 1 ) ;
% evaluate Jacobian in the random vector
H_c = assemble_H( coordinates, elements3, C, R_u, R_v ) ;
J_c = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;


% test problem
% C = randn(2, 1) ;
% H = @(C) [ 2*C(1)*C(2) ; C(2)^2 ; exp(C(1)) ] ;
% J = @(C) [ 2*C(2) , 2*C(1) ; 0 , 2*C(2) ; exp(C(1)) , 0 ] ;
% J_c = J(C) ;
% H_c = H(C) ;

% maximal exponent
E = 16 ; path = zeros(1, E) ;
% store max difference between jacobian and approximation
for e = 1:E
	% approxiamte jacobian
    J_appr = zeros(size(J_c)) ;

    for j = 1:size(J_appr, 2)
    	% perturbation
    	d_C = zeros(size(C)) ;
    	d_C(j) = 10^(-e) ;

    	% H_pdc = H(C + d_C) ;
    	% H_mdc = assemble_H( coordinates, elements3, C - d_C, R_u, R_v ) ;
	    H_pdc = assemble_H( coordinates, elements3, C + d_C, R_u, R_v ) ;
    	for i = 1:size(J_appr, 1)
    		J_appr(i, j) = (H_pdc(i) - H_c(i)) / (10^(-e)) ;
    	end
    end

    path(e) = max(max( J_appr - J_c )) ;
end

% plot convergence
figure
box on
hold on
plot(10.^(-1:-1:-E), path, 'k-', 'linewidth', 1) ;
plot(10.^(-1:-1:-E), 10.^(-1:-1:-E), 'b--', 'linewidth', 1)
hold off
set(gca, 'XScale', 'log', 'YScale', 'log')
legend('Empirical convergence rate', 'Expected convergence rate', 'location', 'best')
title('Convergence of Finite difference approximaiton of the Jacobian')
xlabel('Order of perturbation on C')
ylabel({'Difference between', 'Jacobian and approximation'})