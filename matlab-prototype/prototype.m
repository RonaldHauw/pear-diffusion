%% Group 11 - March 9th 2020
% finite element method for two-dimensional reaction-diffusion system
% based on https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf

clear all

load workspace.mat ;

% compare results with Lammertyn-2003b
% O2    : low  concentration in center, increases towards surface
% CO2   : high concentration in center, decreases towards surface

% LOAD DOMAIN
%
load meshes/HalfCircleMesh.mat
load meshes/HalfCircleMesh_Data.mat
%
% load meshes/HCT_Mesh.mat
% load meshes/HCT_Mesh_Data.mat
%
% load meshes/HCT_Fine_Mesh.mat
% load meshes/HCT_Fine_Mesh_Data.mat
%
coordinates = mesh.Nodes'/30 ;
% coordinates(:, 2) = coordinates(:, 2) + abs(min(coordinates(:, 2))) ;
elements3   = Elements( : , 2:end ) ;
% number of vertices
M           = size(coordinates, 1) ;
% coordinates of mesh vertices
r           = coordinates(:, 1) ;
z           = coordinates(:, 2) ;
% edge information
G1_edges    = InnerBEdges( :, 2:end ) ;
G1_nodes    = InnerBNodes' ;
G2_edges    = OuterBEdges( :, 2:end ) ;
G2_nodes    = OuterBNodes' ;

% INITIAL VALUE
% find solution of linearized system to initialize Newton-Raphson

% K = [ K_u , 0 ; 0 , K_v ]
K = assemble_K( coordinates, elements3, G2_edges, ...
                s_ur, s_vr, s_uz, s_vz, r_u, r_v ) ;
% f = [ f_u ; f_v ]
f = assemble_f( coordinates, G2_edges, ...
                r_u, r_v, C_u_amb, C_v_amb ) ;
% linearization of H = [ H_u(c_u, c_v) ; H_v(c_u, c_v)] around (C_u_amb, C_v_amb)
[J, l] = assemble_H_lin( coordinates, elements3, ...
                         C_u_amb, C_v_amb, R_u, R_v, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;

% C holds the coefficients c_i and c_{M+i}
C = (K+J)\(f-l) ;
% C = K\f ; % check only diffusion


% DOES NOT WORK make first newton step from ambient concentrations
% C_amb = [ C_u_amb*ones(M, 1) ; C_v_amb*ones(M, 1) ] ;
% H_amb = assemble_H( coordinates, elements3, C_amb, R_u, R_v ) ;
% J_amb = assemble_J( coordinates, elements3, C_amb, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
% g_amb = K*C_amb - f + H_amb ; 
% C = C_amb - (K+J_amb)\(g_amb) ;

% C = max(0, C) ;
% C = [ C_u_amb*ones(M, 1) ; C_v_amb*ones(M, 1) ] ;


C(1:M) = C(1:M) + abs(min(C(1:M))) ;
C(M+1:end) = C(M+1:end) + abs(min(C(M+1:end))) ;

% SOLVE NONLINEAR SYSTEM
% Newton-Raphson iteration
for n=1:50
    
    % no need to recompute K and f
    
    % H = [ H_u(c_u, c_v) ; H_v(c_u, c_v)]
    H = assemble_H( coordinates, elements3, ...
                    C, R_u, R_v ) ;
    % Jacobian J = K + dH/dC
    J = assemble_J( coordinates, elements3, ...
                    C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
    
    % Variational
    G = K*C - f + H ;
    
    % Solving one Newton step
    P = (K+J)\G ;
    C = C - P;
    
%     subplot(1, 2, 1)
%     show(elements3,[],coordinates,full( C(1:M) ));
%     title('O_2 concentration (%)')
%     subplot(1, 2, 2)
%     show(elements3,[],coordinates,full( C(M+1:end) ));
%     title('CO_2 concentration (%)')
%     suptitle(join(['Conditions : ', num2str(100*n_u), '% O_2, ', num2str(100*n_v), '% CO_2, and ', num2str(T_cel), '°C' ]))
%     pause(1)
    
    % check for convergence
    norm(P)
    if norm(P) < 10^(-8)
        disp("stop at iteration " + num2str(n) ) ;
        break
    end
end


% GRAPHIC REPRESENTATION OF CONCENTRATIONS

% convert units from mol/m^3 to %
density     = 970 ;     % density of pear in kg/m^3 (see Lammertyn 2003a)
mass_o2     = 3.2e-2 ;  % molar mass of oxygen in kg/mol
mass_co2    = 4.4e-2 ;  % molar mass of carbon dioxide in kg/mol
%
C(1:M)      = 100 * C(1:M)     * mass_o2  / density ;
C(M+1:end)  = 100 * C(M+1:end) * mass_co2 / density ;

figure('position', [300 100 800 500])
subplot(1, 2, 1)
show(elements3,[],coordinates,full( C(1:M) ));
title('O_2 concentration (%)')
subplot(1, 2, 2)
show(elements3,[],coordinates,full( C(M+1:end) ));
title('CO_2 concentration (%)')
suptitle(join(['Conditions : ', num2str(100*n_u), '% O_2, ', num2str(100*n_v), '% CO_2, and ', num2str(T_cel), '°C' ]))