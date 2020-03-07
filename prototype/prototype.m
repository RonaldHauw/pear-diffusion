%% Group 11 - 21st February 2020
% finite element method for two-dimensional nonlinear equation
% based on https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf

clear all

% LOAD DOMAIN
%
 %load mesh/HalfCircleMesh.mat
 %load mesh/HalfCircleMesh_Data.mat
%
load mesh/HCTmesh3.mat
load mesh/HCTmesh3_Data.mat
%
coordinates = Nodes;
r           = coordinates(:, 2) ;
z           = coordinates(:, 3) ;
elements3   = Elements( : , 2:end ) ;
G1_edges    = InnerBEdges( :, 2:end ) ;
G1_nodes    = InnerBNodes' ;
G2_edges    = OuterBEdges( :, 2:end ) ;
G2_nodes    = OuterBNodes' ;


% DECLARE PARAMETERS
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
R_u     = @(c_u, c_v) V_mu * c_u ./ ( (K_mu+c_u) .* (1+c_v/K_mv) ) ;
R_v     = @(c_u, c_v) r_q*R_u(c_u, c_v) +  V_mfv ./ (1+c_u/K_mfu) ;
% derivative of respiration
dR_u_u  = @(c_u, c_v) V_mu ./ (K_mu+c_u) ./ (1+c_v/K_mv) .* (1-c_u./(K_mu+c_u)) ;
dR_u_v  = @(c_u, c_v) - 1/K_mv * V_mu * c_u / (K_mu+c_u) / (1+c_v/K_mv)^2 ;
dR_v_u  = @(c_u, c_v) r_q*dR_u_u(c_u, c_v) - 1/K_mfu * V_mfv / (1+c_u/K_mfu)^2 ;
dR_v_v  = @(c_u, c_v) r_q*dR_u_v(c_u, c_v) ;
% linearization of respiration around (c_u_amb, C_v_amb)
R_u_lin = @(c_u, c_v) R_u(C_u_amb, C_v_amb) + dR_u_u(C_u_amb, C_v_amb)*(c_u-C_u_amb) ...
                                            + dR_u_v(C_u_amb, C_v_amb)*(c_v-C_v_amb) ;
R_v_lin = @(c_u, c_v) R_v(C_u_amb, C_v_amb) + dR_v_u(C_u_amb, C_v_amb)*(c_u-C_u_amb) ...
                                            + dR_v_v(C_u_amb, C_v_amb)*(c_v-C_v_amb) ;

% Nodes where no dirichlet is imposed (value is "free")
% FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));

% INITIAL VALUE
% C holds the coefficients c_i and c_{M+i}
C = ones(2*M,1);


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

% set up linear system to solve
A = K + H ;
b = f - l ;
%C = A\b ;


% Newton-Raphson iteration
for n=1:200
    
    % K = [ K_u , 0 ; 0 , K_v ]
    K = assemble_K( coordinates, elements3, G2_edges, ...
                    s_ur, s_vr, s_uz, s_vz, r_u, r_v );
    % f = [ f_u ; f_v ]
    f = assemble_f( coordinates, G2_edges, ...
                    r_u, r_v, C_u_amb, C_v_amb ) ;
    % H = [ H_u(c_u, c_v) ; H_v(c_u, c_v)]
    H = assemble_H( coordinates, elements3, ...
                    C, R_u, R_v ) ;
    % Jacobian J = K + dH/dC
    J = assemble_J( coordinates, elements3, ...
                    C, K, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
    
    % Variational
    G = K*C - f + H ;
    
    % Solving one Newton step
    P = J\G ;
    C = C - 0.05*P;
    
    % check for convergence
    norm(P)
    if norm(P) < 10^(-8)
        disp("stop at iteration " + num2str(n) ) ;
        break
    end
end

% graphic representation
figure()
subplot(1, 2, 1)
show(elements3,[],coordinates(:, 2:3),C(1:M)); % full(C(1:M))
title('oxygen')
subplot(1, 2, 2)
show(elements3,[],coordinates(:, 2:3),C(M+1:end));
title('carbon dioxide')

%%

% Assemble diffusion matrix K = [ K_u , 0 ; 0 , K_v ]
function K = assemble_K( coordinates, elements3, G2_edges, s_ur, s_vr, s_uz, s_vz, r_u, r_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements
    % G2_edges          index of vertices that form the edges ofouter boundary
    
    % extract useful variables
    M = size(coordinates, 1) ;
    
    % CHANGED: commented
    %r = coordinates(:, 2) ;
    %z = coordinates(:, 3) ;

    K = zeros( 2*M, 2*M );
    for t = elements3'
        % area of element (can be positive or negative)
        omega = det([ ones(1,3) ; coordinates(t, 2:3)' ]) / 2;
        
        % CHANGED: added
        r = coordinates(t, 2)'; 
        z = coordinates(t, 3)';


        % sum of r-coordinates
        sum_r = sum(coordinates(t, 2), 1); 
        
        
        % for j different from i
        C_12 = 1/6 * 1/2/omega * [ (z(1)-z(3))*(z(3)-z(2)) ; (r(1)-r(3))*(r(3)-r(2))] ;
        C_23 = 1/6 * 1/2/omega * [ (z(2)-z(1))*(z(1)-z(3)) ; (r(2)-r(1))*(r(1)-r(3))] ;
        C_13 = 1/6 * 1/2/omega * [ (z(1)-z(2))*(z(2)-z(3)) ; (r(1)-r(2))*(r(2)-r(3))] ;
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
        C_11 = 1/6 * 1/2/omega * [ (z(2)-z(3))^2 ; (r(2)-r(3))^2] ;
        C_22 = 1/6 * 1/2/omega * [ (z(1)-z(3))^2 ; (r(1)-r(3))^2] ;
        C_33 = 1/6 * 1/2/omega * [ (z(1)-z(2))^2 ; (r(1)-r(2))^2] ;
        
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
    
    % CHANGED: added
    % bug: waarom enkel r component hier? 
    r = coordinates(:, 2) ;
    z = coordinates(:, 3) ;
    for e = G2_edges'
        % length of edge
        len = norm(diff(coordinates(e, 2:3)), 2) ;
        
        
        % compute two different terms
        parallel_term_1     = 1./12 * len * ( 3*r(e(1)) +   r(e(2)) ) ;
        parallel_term_2     = 1./12 * len * (   r(e(1)) + 3*r(e(2)) ) ;
        cross_term          = 1./12 * len * (   r(e(1)) +   r(e(2)) ) ;
        
        % in K_u
        K( e(1),   e(1) )   = K( e(1), e(1) )     + r_u * parallel_term_1 ;
        K( e(1),   e(2) )   = K( e(1), e(2) )     + r_u * cross_term ;
        K( e(2),   e(1) )   = K( e(2), e(1) )     + r_u * cross_term ; 
        K( e(2),   e(2) )   = K( e(2), e(2) )     + r_u * parallel_term_2 ;
        % in K_v
        K( M+e(1), M+e(1) ) = K( M+e(1), M+e(1) ) + r_v * parallel_term_1 ;
        K( M+e(1), M+e(2) ) = K( M+e(1), M+e(2) ) + r_v * cross_term ;
        K( M+e(2), M+e(1) ) = K( M+e(2), M+e(1) ) + r_v * cross_term ; 
        K( M+e(2), M+e(2) ) = K( M+e(2), M+e(2) ) + r_v * parallel_term_2 ;
    end

end


% Assemble f = [ f_u ; f_v ]
function f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb )
    % coordinates       coordinates of vertices of mesh
    % G2_edges          index of vertices that form the edges ofouter boundary 

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 2) ;
    z = coordinates(:, 3) ;
    
    f = zeros( 2*M, 1 ) ;
    for e = G2_edges'
        % length of edge
        len = norm(diff(coordinates(e, :)), 2) ;
        
        % compute two terms
        term_1 = 1/6 * len * ( 2*r(e(1)) +   r(e(2)) ) ;
        term_2 = 1/6 * len * (   r(e(1)) + 2*r(e(2)) ) ;
        
        % in f_u
        f( e(1) )   = f( e(1) )   + r_u * C_u_amb * term_1 ;
        f( e(2) )   = f( e(2) )   + r_u * C_u_amb * term_2 ;
        % in f_v
        f( M+e(1) ) = f( M+e(1) ) + r_v * C_v_amb * term_1 ;
        f( M+e(2) ) = f( M+e(2) ) + r_v * C_v_amb * term_2 ;
    end
end


% Assemble H = [ H_u ; H_v]
function H = assemble_H( coordinates, elements3, C, R_u, R_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 2) ;
    
    H = zeros( 2*M, 1 ) ;
    for t = elements3'
        % area of element
        area = abs(det([ ones(1,3) ; coordinates(t, 2:3)' ])) / 2 ;
        
        % in H_u
        H( t )   = H( t )   + 1/3*area * r(t) .* R_u(C(t), C(M+t)) ;
        % in H_v
        H( M+t ) = H( M+t ) - 1/3*area * r(t) .* R_v(C(t), C(M+t)) ;
        
    end
    H_= H
end


function [H, l] = assemble_H_lin( coordinates, elements3, C_u_amb, C_v_amb, R_u, R_v, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements
    
    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 2) ;
    
    H = zeros(2*M, 2*M) ;
    l = zeros(2*M, 1) ;
    for t = elements3'
        % area of element
        area = abs(det([ ones(1,3) ; coordinates(t, 2:3)' ])) / 2 ;
        
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
    r = coordinates(:, 2) ;

    J = zeros( 2*M, 2*M ) ;
    for i = 1:M
        % sum of areas of elements arount vertex i
        T = elements3(any(elements3==i, 2), :) ;
        s = 0 ;
        
        for t = T'
            s = s + abs(det([ ones(1,3) ; coordinates(t, 2:3)' ])) / 2 ;
            
        end
  
        
        
        
        % part derivative of H_u to C
        J( i, i )     =  s/3. * r(i) * dR_u_u(C(i), C(M+i)) ;
        J( i, M+i )   =  s/3. * r(i) * dR_u_v(C(i), C(M+i)) ;
                
        % part derivative of H_v to C
        J( M+i, i )   = -s/3 * r(i) * dR_v_u(C(i), C(M+i)) ;
        J( M+i, M+i ) = -s/3 * r(i) * dR_v_v(C(i), C(M+i)) ;
        

        
    end
    J = J + K ;
end
