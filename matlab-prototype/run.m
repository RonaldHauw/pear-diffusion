%% Group 11 - March 9th 2020
% finite element method for two-dimensional reaction-diffusion system
% based on https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf

clear all
clc

% choose 
T_cel   = 25 ;      % degrees in celcius	
n_u     = 0.208 ;   % concentration O2  percentage in 0 < . < 1
n_v     = 0.0 ;       % concentration CO2 percentage in 0 < . < 1
%
workspace ;
load workspace.mat ;

% LOAD DOMAIN
%
% mfilepath=fileparts(which(mfilename));
% addpath(fullfile(mfilepath,'../1-pear-diffusion'));
addpath( '../grids' )
% load HalfCircleMesh.mat
% load HalfCircleMesh_Data.mat
%
load HCT_Mesh.mat
load HCT_Mesh_Data.mat
%
% load HCT_Fine_Mesh.mat
% load HCT_Fine_Mesh_Data.mat
%
coordinates = mesh.Nodes'/50 ;
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

% parameters of homotopy continuation
dt = 0.05 ;
maxit = 50 ;

% intiialize concentrations
C = zeros(2*M, 1) ;
% K = [ K_u , 0 ; 0 , K_v ]
K = assemble_K( coordinates, elements3, G2_edges, s_ur, s_vr, s_uz, s_vz, r_u, r_v ) ;
% f = [ f_u ; f_v ]
f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb ) ;

% perform hopotopy continuation
for t = 0:dt:1
    
    % prediction step
    % nonlinearity H = [ H_u(C) ; H_v(C) ]
    H       = assemble_H( coordinates, elements3, C, R_u, R_v ) ;
    % Jacobian J = dH/dC
    J       = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
    % update concentrations with forward euler
    dC_dt   = ( K + t*J ) \ ( -H ) ;
    C       = C + dt * dC_dt ;
    
    
    % correction step with Newton method
    for n=1:maxit

        % no need to recompute K and f

        % nonlinearity H = [ H_u(C) ; H_v(C) ]
        H = assemble_H( coordinates, elements3, C, R_u, R_v ) ;
        % Jacobian J = dH/dC
        J = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;

        % Variational
        G = K*C - f + t*H ;

        % solving one Newton step (J_G)^-1 * G
        P = ( K + t*J ) \ G ;
        
        % check for convergence
        if norm(P) < 10^(-12)
            disp("t = " + num2str(t) + "    stop at itr " + num2str(n) + "    with resid " + num2str(norm(P)) ) ;
            break
        end
        
        % backtracking
        b = 1 ;
        for k = 1:50
            % store proposed new concentrations
            temp = C - b*P ;
            
            % recompute Variational
            H = assemble_H( coordinates, elements3, temp, R_u, R_v ) ;
            res = K*temp - f + t*H ;
            
            % check if new concentrations indeed reduce the Variational
            if ( norm(res) > norm(G) )
                b = b/2 ;
            else
                break
            end
        end
        
        % update estimate for concentrations
        C = C - b*P;
    end
end


% GRAPHIC REPRESENTATION OF CONCENTRATIONS

% % convert units from mol/m^3 to %
% density     = 970 ;     % density of pear in kg/m^3 (see Lammertyn 2003a)
% mass_o2     = 3.2e-2 ;  % molar mass of oxygen in kg/mol
% mass_co2    = 4.4e-2 ;  % molar mass of carbon dioxide in kg/mol
% %
% C(1:M)      = 100 * C(1:M)     * mass_o2  / density ;
% C(M+1:end)  = 100 * C(M+1:end) * mass_co2 / density ;
%%
figure('position', [300 100 450 400])
subplot(1, 2, 1)
show(elements3,[],coordinates,full( C(1:M) ), [0, 20], C_u_amb);
title('O_2 concentration (%)', 'FontSize', 12)
subplot(1, 2, 2)
show(elements3,[],coordinates,full( C(M+1:end) ), [0, 6], C_v_amb);
title('CO_2 concentration (%)', 'FontSize', 12)
% suptitle( join(['Conditions : ', num2str(100*n_u), '% O_2, ', num2str(100*n_v), '% CO_2, and ', num2str(T_cel), '°C' ]), 'FontSize', 12 )