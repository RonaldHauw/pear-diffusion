%% Group 11 - 21st February 2020
% finite element method for two-dimensional nonlinear equation
% based on https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf

clear all



% LOAD DOMAIN
%
load mesh/HalfCircleMesh.mat
load mesh/HalfCircleMesh_Data.mat
%
% load mesh/HCT_Mesh.mat
% load mesh/HCT_Mesh_Data.mat
%
% load mesh/HCT_Fine_Mesh.mat
% load mesh/HCT_Fine_Mesh_Data.mat
%
coordinates = mesh.Nodes' ;
r           = coordinates(:, 1) ;
z           = coordinates(:, 2) ;
elements3   = Elements( : , 2:end ) ;
G1_edges    = InnerBEdges( :, 2:end ) ;
G1_nodes    = InnerBNodes' ;
G2_edges    = OuterBEdges( :, 2:end ) ;
G2_nodes    = OuterBNodes' ;



% DECLARE PARAMETERS and functions
% variable
T_cel           = 20 ;      % degrees in celcius	
n_u             = 0.208 ;   % percentage in 0 < . < 1
n_v             = 0 ;       % percentage in 0 < . < 1
% fixed
M               = size(coordinates, 1) ;
%
s_ur            = 2.8e-10 ;
s_uz            = 1.1e-9 ;
s_vr            = 2.32e-9 ;
s_vz            = 6.97e-9 ;
%
T_ref           = 293.15 ;
R_g             = 8.314 ;
V_mu_ref        = 2.39e-4 ;
E_a_vmu_ref     = 80200 ;
T               = T_cel + 273.15 ;
V_mu            = V_mu_ref * exp( E_a_vmu_ref/R_g * ( 1/T_ref-1/T ) ) ;
%
V_mfv_ref       = 1.61e-4 ;
E_a_vmfv_ref    = 56700 ;
V_mfv           = V_mfv_ref * exp( E_a_vmfv_ref/R_g * ( 1/T_ref-1/T ) ) ;
%
K_mu            = 0.4103 ;
K_mv            = 27.2438 ;
K_mfu           = 0.1149 ;
%
r_q             = 0.97 ;
r_u             = 7e-7 ; 
r_v             = 7.5e-7 ;
%
p_atm           = 101300 ;
C_u_amb         = p_atm*n_u/R_g/T ;
C_v_amb         = p_atm*n_v/R_g/T ;
% respiration
R_u     		= @(c_u, c_v) K_mv * V_mu * c_u ./ ( (K_mu+c_u) .* (K_mv+c_v) ) ;
R_v     		= @(c_u, c_v) r_q*R_u(c_u, c_v) +  V_mfv ./ (1+c_u/K_mfu) ;
% derivative of respiration
dR_u_u  		= @(c_u, c_v) K_mv * V_mu ./ (K_mu+c_u) ./ (K_mv+c_v) .* (1-c_u./(K_mu+c_u)) ;
dR_u_v  		= @(c_u, c_v) - K_mv * V_mu * c_u ./ (K_mu+c_u) ./ (K_mv+c_v).^2 ;
dR_v_u  		= @(c_u, c_v) r_q*dR_u_u(c_u, c_v) - K_mfu * V_mfv ./ (K_mfu+c_u).^2 ;
dR_v_v  		= @(c_u, c_v) r_q*dR_u_v(c_u, c_v) ;
% linearization of respiration around (c_u_amb, C_v_amb)
R_u_lin 		= @(c_u, c_v) R_u(C_u_amb, C_v_amb) + dR_u_u(C_u_amb, C_v_amb)*(c_u-C_u_amb) ...
                                           		    + dR_u_v(C_u_amb, C_v_amb)*(c_v-C_v_amb) ;
R_v_lin 		= @(c_u, c_v) R_v(C_u_amb, C_v_amb) + dR_v_u(C_u_amb, C_v_amb)*(c_u-C_u_amb) ...
                                          		    + dR_v_v(C_u_amb, C_v_amb)*(c_v-C_v_amb) ;



% INITIAL VALUE
% find solution of linearized system to initialize Newton-Raphson

% K = [ K_u , 0 ; 0 , K_v ]
K = assemble_K( coordinates, elements3, G2_edges, ...
                s_ur, s_vr, s_uz, s_vz, r_u, r_v ) ;
% f = [ f_u ; f_v ]
f = assemble_f( coordinates, G2_edges, ...
                r_u, r_v, C_u_amb, C_v_amb ) ;
% linearization of H = [ H_u(c_u, c_v) ; H_v(c_u, c_v)] around (C_u_amb, C_v_amb)
[H, l] = assemble_H_lin( coordinates, elements3, ...
                         C_u_amb, C_v_amb, R_u, R_v, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
% minus sign in front of second linearization in R_v is already accounted for

% C holds the coefficients c_i and c_{M+i}
C = (K+H)\(f-l) ;

% C = [ C_u_amb*ones(M, 1) ; C_v_amb*ones(M, 1) ] ;


% SOLVE NONLINEAR SYSTEM
% Newton-Raphson iteration
% for n=1:50
%     
%     % K = [ K_u , 0 ; 0 , K_v ]
%     K = assemble_K( coordinates, elements3, G2_edges, ...
%                     s_ur, s_vr, s_uz, s_vz, r_u, r_v ) ;
%     % f = [ f_u ; f_v ]
%     f = assemble_f( coordinates, G2_edges, ...
%                     r_u, r_v, C_u_amb, C_v_amb ) ;
%     % H = [ H_u(c_u, c_v) ; H_v(c_u, c_v)]
%     H = assemble_H( coordinates, elements3, ...
%                     C, R_u, R_v ) ;
%     % Jacobian J = K + dH/dC
%     J = assemble_J( coordinates, elements3, ...
%                     C, K, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
%     
%     % Variational
%     G = K*C - f + H ;
%     
%     % Solving one Newton step
%     P = J\G ;
%     C = C - P;
%     
%     % check for convergence
%     norm(P)
%     if norm(P) < 10^(-10)
%         disp("stop at iteration " + num2str(n) ) ;
%         break
%     end
% end


% GRAPHIC REPRESENTATION OF CONCENTRATIONS

% convert units from mol/m^3 to %
density     = 970 ;     % density of pear in kg/m^3 (see Lammertyn 2003a)
mass_o2     = 3.2e-2 ;  % molar mass of oxigen in kg/mol
mass_co2    = 4.4e-2 ;  % molar mass of carbon dioxide in kg/mol
% convert
C(1:M)     = C(1:M)     * mass_o2  / density ;
C(M+1:end) = C(M+1:end) * mass_co2 / density ;

figure('position', [300 100 800 500])
subplot(1, 2, 1)
show(elements3,[],coordinates,full(C(1:M)));
title('O_2 concentration (%)')
subplot(1, 2, 2)
show(elements3,[],coordinates,full(C(M+1:end)));
title('CO_2 concentration (%)')

suptitle(join(['Conditions : ', num2str(100*n_u), '% O_2, ', num2str(100*n_v), '% CO_2, and ', num2str(T_cel), '°C' ]))


% CHECK FOR CORRECTNESS OF RESULTS
% norm(C(unique(G2_edges))   - C_u_amb*mass_o2/density*ones(size(unique(G2_edges), 1), 1))
% norm(C(unique(M+G2_edges)) - C_v_amb*mass_co2/density*ones(size(unique(G2_edges), 1), 1))


% FUNCTIONS
% Assemble diffusion matrix K = [ K_u , 0 ; 0 , K_v ]
function K = assemble_K( coordinates, elements3, G2_edges, s_ur, s_vr, s_uz, s_vz, r_u, r_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements
    % G2_edges          index of vertices that form the edges ofouter boundary
    
    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;
    z = coordinates(:, 2) ;

    K = zeros( 2*M, 2*M );
    for t = elements3'
                
        % area of element (can be positive or negative)
        omega = det([ ones(1,3) ; coordinates(t, :)' ]) / 2 ;
		% omega_12 = det([ ones(1,3) ; coordinates(t, :)' ]) / 2
		% omega_21 = det([ ones(1,3) ; coordinates(t([1, 3, 2]), :)' ]) / 2

        % sum of r-coordinates
        sum_r = sum(coordinates(t, 1), 1) ;
        
        % coordinates of elements vertices
        
        
        % for j different from i
        C_12 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(3)))*(z(t(3))-z(t(2))) ; ...
                                   (r(t(1))-r(t(3)))*(r(t(3))-r(t(2)))] ;
                               
        C_23 = 1/6 * 1/2/omega * [ (z(t(2))-z(t(1)))*(z(t(1))-z(t(3))) ; ...
                                   (r(t(2))-r(t(1)))*(r(t(1))-r(t(3)))] ;
                               
        C_13 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(2)))*(z(t(2))-z(t(3))) ; ...
                                   (r(t(1))-r(t(2)))*(r(t(2))-r(t(3)))] ;
        %
        K(t(1),   t(2))   = K(t(1),   t(2))   + [s_ur, s_uz] * C_12 * sum_r ;
        K(t(2),   t(1))   = K(t(2),   t(1))   + [s_ur, s_uz] * C_12 * sum_r ;
        %
        K(t(2),   t(3))   = K(t(2),   t(3))   + [s_ur, s_uz] * C_23 * sum_r ;
        K(t(3),   t(2))   = K(t(3),   t(2))   + [s_ur, s_uz] * C_23 * sum_r ;
        %
        K(t(1),   t(3))   = K(t(1),   t(3))   + [s_ur, s_uz] * C_13 * sum_r ;
        K(t(3),   t(1))   = K(t(3),   t(1))   + [s_ur, s_uz] * C_13 * sum_r ;

        K(M+t(1), M+t(2)) = K(M+t(1), M+t(2)) + [s_vr, s_vz] * C_12 * sum_r ;
        K(M+t(2), M+t(1)) = K(M+t(2), M+t(1)) + [s_vr, s_vz] * C_12 * sum_r ;
        %
        K(M+t(2), M+t(3)) = K(M+t(2), M+t(3)) + [s_vr, s_vz] * C_23 * sum_r ;
        K(M+t(3), M+t(2)) = K(M+t(3), M+t(2)) + [s_vr, s_vz] * C_23 * sum_r ;
        %
        K(M+t(1), M+t(3)) = K(M+t(1), M+t(3)) + [s_vr, s_vz] * C_13 * sum_r ;
        K(M+t(3), M+t(1)) = K(M+t(3), M+t(1)) + [s_vr, s_vz] * C_13 * sum_r ;
        
        % for j equal i
        C_11 = 1/6 * 1/2/omega * [ (z(t(2))-z(t(3)))^2 ; (r(t(2))-r(t(3)))^2] ;
        C_22 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(3)))^2 ; (r(t(1))-r(t(3)))^2] ;
        C_33 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(2)))^2 ; (r(t(1))-r(t(2)))^2] ;
        %
        K(t(1),   t(1))   = K(t(1),   t(1))   + [s_ur, s_uz] * C_11 * sum_r ;
        K(t(2),   t(2))   = K(t(2),   t(2))   + [s_ur, s_uz] * C_22 * sum_r ;
        K(t(3),   t(3))   = K(t(3),   t(3))   + [s_ur, s_uz] * C_33 * sum_r ;
        %
        K(M+t(1), M+t(1)) = K(M+t(1), M+t(1)) + [s_vr, s_vz] * C_11 * sum_r ;
        K(M+t(2), M+t(2)) = K(M+t(2), M+t(2)) + [s_vr, s_vz] * C_22 * sum_r ;
        K(M+t(3), M+t(3)) = K(M+t(3), M+t(3)) + [s_vr, s_vz] * C_33 * sum_r ;

    end
    % add terms for vertices on outer boundary
    % assume outer boundary goes from bottom to top
    for e = G2_edges'
        % length of edge
        len = norm(diff(coordinates(e, :)), 2) ;
        
        % compute two different terms
        parallel_term_1     = 1/12 * len * ( 3*r(e(1)) +   r(e(2)) ) ;
        parallel_term_2     = 1/12 * len * (   r(e(1)) + 3*r(e(2)) ) ;
        cross_term          = 1/12 * len * (   r(e(1)) +   r(e(2)) ) ;
        
        % in K_u
        % K( e(1),   e(1) )   = K( e(1), e(1) )     + r_u * parallel_term_1 ;
        % K( e(1),   e(2) )   = K( e(1), e(2) )     + r_u * cross_term ;
        % K( e(2),   e(1) )   = K( e(2), e(1) )     + r_u * cross_term ; 
        % K( e(2),   e(2) )   = K( e(2), e(2) )     + r_u * parallel_term_2 ;
        K( e,   e )   = K( e,   e )   + [ r_u * parallel_term_1 , r_u * cross_term ; ...
        							      r_u * cross_term      , r_u * parallel_term_2 ] ;
        % in K_v
        % K( M+e(1), M+e(1) ) = K( M+e(1), M+e(1) ) + r_v * parallel_term_1 ;
        % K( M+e(1), M+e(2) ) = K( M+e(1), M+e(2) ) + r_v * cross_term ;
        % K( M+e(2), M+e(1) ) = K( M+e(2), M+e(1) ) + r_v * cross_term ; 
        % K( M+e(2), M+e(2) ) = K( M+e(2), M+e(2) ) + r_v * parallel_term_2 ;
        K( M+e, M+e ) = K( M+e, M+e ) + [ r_v * parallel_term_1 , r_v * cross_term ; ...
        							      r_v * cross_term      , r_v * parallel_term_2 ] ;
    end

end


% Assemble f = [ f_u ; f_v ]
function f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb )
    % coordinates       coordinates of vertices of mesh
    % G2_edges          index of vertices that form the edges ofouter boundary 

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;
    
    f = zeros( 2*M, 1 ) ;
    for e = G2_edges'
        % length of edge
        len = norm(diff(coordinates(e, :)), 2) ;
        
        % compute two terms
        term_1 = 1/6 * len * ( 2*r(e(1)) +   r(e(2)) ) ;
        term_2 = 1/6 * len * (   r(e(1)) + 2*r(e(2)) ) ;
        
        % in f_u
        %f( e(1) )   = f( e(1) )   + r_u * C_u_amb * term_1 ;
        %f( e(2) )   = f( e(2) )   + r_u * C_u_amb * term_2 ;
        f( e )   = f( e )   + r_u * C_u_amb * [ term_1 ; term_2 ] ;
        % in f_v
        %f( M+e(1) ) = f( M+e(1) ) + r_v * C_v_amb * term_1 ;
        %f( M+e(2) ) = f( M+e(2) ) + r_v * C_v_amb * term_2 ;
        f( M+e ) = f( M+e ) + r_v * C_v_amb * [ term_1 ; term_2 ] ;
    end
end


% Assemble H = [ H_u ; H_v]
function H = assemble_H( coordinates, elements3, C, R_u, R_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;
    
    H = zeros( 2*M, 1 ) ;
    for t = elements3'
        % area of element
        area = abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        
        % in H_u
        H( t )   = H( t )   + 1/3 * area * r(t) .* R_u(C(t), C(M+t)) ;
        % in H_v
        H( M+t ) = H( M+t ) - 1/3 * area * r(t) .* R_v(C(t), C(M+t)) ;
        
    end
end


function [H, l] = assemble_H_lin( coordinates, elements3, C_u_amb, C_v_amb, R_u, R_v, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements
    
    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;
    
    H = zeros(2*M, 2*M) ;
    l = zeros(2*M, 1) ;
    for t = elements3'
        % area of element
        area = abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        
        % contribution of linearized R_u
        l(t)   = l(t)   + 1/3 * area * r(t) .* ( R_u(C_u_amb, C_v_amb) ...
                                               - C_u_amb*dR_u_u(C_u_amb, C_v_amb) ...
                                               - C_v_amb*dR_u_v(C_u_amb, C_v_amb)) ;
        % contribution of linearized R_v
        l(M+t) = l(M+t) - 1/3 * area * r(t) .* ( R_v(C_u_amb, C_v_amb) ...
                                               - C_u_amb*dR_v_u(C_u_amb, C_v_amb) ...
                                               - C_v_amb*dR_v_v(C_u_amb, C_v_amb)) ;

        for i = t'
            % contribution from R_u
            H(i, i)     = H(i, i)     + 1/3 * area * r(i) .* dR_u_u(C_u_amb, C_v_amb) ;
            H(i, M+i)   = H(i, M+i)   + 1/3 * area * r(i) .* dR_u_v(C_u_amb, C_v_amb) ;
            % contribution from R_v
            H(M+i, i)   = H(M+i, i)   - 1/3 * area * r(i) .* dR_v_u(C_u_amb, C_v_amb) ; 
            H(M+i, M+i) = H(M+i, M+i) - 1/3 * area * r(i) .* dR_v_v(C_u_amb, C_v_amb) ;
        end
        
    end
    
end


% Assemble Jacobian J = K + dH/dC
function J = assemble_J( coordinates, elements3, C, K, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;

    J = zeros( 2*M, 2*M ) ;
    for i = 1:M
        % sum of areas of elements arount vertex i
        T = elements3(any(elements3==i, 2), :) ;
        s = 0 ;
        for t = T'
            s = s + abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        end
        
        % part derivative of H_u to C
        J( i, i )     =  s/3 * r(i) * dR_u_u(C(i), C(M+i)) ;
        J( i, M+i )   =  s/3 * r(i) * dR_u_v(C(i), C(M+i)) ;
        
        % part derivative of H_v to C
        J( M+i, i )   = -s/3 * r(i) * dR_v_u(C(i), C(M+i)) ;
        J( M+i, M+i ) = -s/3 * r(i) * dR_v_v(C(i), C(M+i)) ;
        
    end
    J = J + K ;
end