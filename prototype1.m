%% Group 11 - 17th February 2020
% quick and dirty finite element method for two-dimensional nonlinear equation
% based on https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf

% clear all


% LOAD DOMAIN
% Initialisation
% load coordinates.dat; coordinates(:,1)=[];
% load elements3.dat; elements3(:,1)=[];
% eval('load neumann.dat; neumann(:,1) = [];','neumann=[];');
% load dirichlet.dat; dirichlet(:,1) = [];
load HalfCircleMesh.mat
load HalfCircleMesh_Data.mat
coordinates = mesh.Nodes' ;
r           = coordinates(:, 1) ;
z           = coordinates(:, 2) ;
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



% Nodes where no dirichlet is imposed (value is "free")
% FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));

% INITIAL VALUE
% C holds the coefficients c_i and c_{M+i}
C = ones(2*size(coordinates,1),1);
% U(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));


% Newton-Raphson iteration
for n=1:50
  
%     % Assembly of DJ(U)
%     A = sparse(size(coordinates,1),size(coordinates,1));
%   % for each element (triangle in elements3)
%   for j = 1:size(elements3,1)
%     A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
% 	+ localdj( coordinates(elements3(j,:),:) , U(elements3(j,:)) ) ;
%   end


    % Assemble K = [ K_u 0 ; 0 K_v ]
    K = sparse( 2*M, 2*M );
    % FOR BETTER EFFICIENCY
    %   - check if can assemble K by elements instead of by vertices
    %   - check if K is symmetric to avoid computing twice the same thing
    % for each vertex (there are M vertices)
    for i = 1:M
        % all vertices j that are connected to vertex i
        J = unique(elements3(any(elements3==i, 2), :)) ;
        for j = J(:)'
            % all elements that have i and j as vertices
            T = elements3(and(any(elements3==i, 2), any(elements3==j, 2)), :) ;
            % compute contribution of each element
            for t = T'
                % area of element t
                area = det([ ones(1,3) ; coordinates(t, :)' ]) / 2 ;
                if j == i
                    % remove node i from element t
                    t(t==i) = [] ;
                    % compute two fraction terms
                    fz = diff(z(t))^2/2/area ;
                    fr = diff(r(t))^2/2/area ;
                else
                    % find vertex of element t that is neither i nor j
                    A = t(and(t~=i, t~=j)) ;
                    % compute two fraction terms
                    fz = diff(z([A, i]))*diff(z([j, A]))/2/area ;
                    fr = diff(r([A, j]))*diff(r([i, A]))/2/area ;                  
                end
                % compute coefficient C_jt
                C_jt_u = s_ur*fz + s_uz*fr ;
                C_jt_v = s_vr*fz + s_vz*fr ;
                % update matrix entry with contribution of element t
                K(i, j)     = K(i, j)     + C_jt_u * sum(r(t)) ;
                K(M+i, M+j) = K(M+i, M+j) + C_jt_v * sum(r(t)) ;
            end                
        end
    end
    % add terms for vertices on outer boundary
    % assume outer boundary goes from bottom to top
    for e = G2_edges'
        % compute two different terms
        parallel_term       = ( r(e(2))-r(e(1)) )*( 3*r(e(2))+r(e(1)) ) /12 ;
        cross_term          = ( r(e(2))^2 - r(e(1))^2 ) /12 ;
        
        % in K_u
        K( e(1),   e(1) )   = K( e(1), e(1) ) + r_u * parallel_term ;
        K( e(1),   e(2) )   = K( e(1), e(2) ) + r_u * cross_term ;
        K( e(2),   e(1) )   = K( e(2), e(1) ) + r_u * cross_term ; 
        K( e(2),   e(2) )   = K( e(2), e(2) ) + r_u * parallel_term ;
        % in K_v
        K( M+e(1), M+e(1) ) = K( M+e(1), M+e(1) ) + r_v * parallel_term ;
        K( M+e(1), M+e(2) ) = K( M+e(1), M+e(2) ) + r_v * cross_term ;
        K( M+e(2), M+e(1) ) = K( M+e(2), M+e(1) ) + r_v * cross_term ; 
        K( M+e(2), M+e(2) ) = K( M+e(2), M+e(2) ) + r_v * parallel_term ;
    end
     
    
    % Assemble f = [ f_u ; f_v ]
    f = zeros( 2*M, 1 ) ;
    for e = G2_edges'
        % compute two terms
        parallel_1  = ( r(e(2))-r(e(1)) ) * ( r(e(2))+2*r(e(1)) ) /6 ;
        parallel_2  = ( r(e(2))-r(e(1)) ) * ( r(e(1))+2*r(e(2)) ) /6 ;
        
        % in f_u
        f( e(1) )   = f( e(1) )   + r_u*C_u_amb * parallel_1 ;
        f( e(2) )   = f( e(2) )   + r_u*C_u_amb * parallel_2 ;
        % in f_v
        f( M+e(1) ) = f( M+e(1) ) + r_v*C_v_amb * parallel_1 ;
        f( M+e(2) ) = f( M+e(2) ) + r_v*C_v_amb * parallel_2 ;
    end
    
    
    % Assemble H = [ H_u ; H_v]
    H = zeros( 2*M, 1 ) ;
    for t = elements3'
        area = det([ ones(1,3) ; coordinates(t, :)' ]) / 2 ;
        
        % in H_u
        H( t )   = H( t )   + 1/3*area * V_mu * r(t) .* C(t) ./ ( (K_mu+C(t)) .* (1+C(M+t)/K_mv) ) ;
        % in H_v
        H( M+t ) = H( M+t ) + 1/3*area * r(t) .* ( r_q*V_mu*C(t) ./ ( (K_mu+C(t)) .* (1+C(M+t)/K_mv) ) ...
                            + V_mfv ./ (1+C(t)/K_mfu) );
    end
    
    % Variational
    G = K*C - f + H ;
    
    
    % Assemble Jacobian J = K + dH/dC
    J = zeros( 2*M, 2*M ) ;
    for i = 1:M
        % sum of areas of elements arount vertex i
        T = elements3(and(any(elements3==i, 2), any(elements3==j, 2)), :) ;
        s = 0 ;
        for t = T'
            s = s + det([ ones(1,3) ; coordinates(t, :)' ]) / 2 ;
        end
        
        % part derivative of H_u to C
        J( i, i )     =   s/3*r(i)*V_mu/(K_mu+C(i))/(1+C(M+i)/K_mv) * (1-C(i)/(K_mu+C(i))) ;
        J( i, M+i )   = - s/3/K_mv * r(i)*V_mu*C(i) / (K_mu+C(i)) / (1+C(M+i)/K_mv)^2 ;
        
        % part derivative of H_v to C
        J( M+i, i )   =   r_q * J( i, i ) - s/3 * V_mfv / K_mfu / (1+C(i)/K_mfu)^2 ;
        J( M+i, M+i)  =   r_q * J( i, M+i ) ;
        
    end
    J = J + K ;
    
    
%   % Assembly of J(U)
%   b = sparse(size(coordinates,1),1);
%   for j = 1:size(elements3,1)
%     b(elements3(j,:)) = b(elements3(j,:)) ...
%    	+ localj(coordinates(elements3(j,:),:),U(elements3(j,:)));
%   end
%   
%   % Volume Forces
%   for j = 1:size(elements3,1)
%     b(elements3(j,:)) = b(elements3(j,:)) + ...
% 	det([1 1 1; coordinates(elements3(j,:),:)']) * ...
% 	f(sum(coordinates(elements3(j,:),:))/3)/6;
%   end
%   
%   % Neumann conditions
%   for j = 1 : size(neumann,1)
%     b(neumann(j,:))=b(neumann(j,:)) - norm(coordinates(neumann(j,1),:)- ...
% 	  coordinates(neumann(j,2),:))*g(sum(coordinates(neumann(j,:),:))/2)/2;
%   end
%   
%   % Dirichlet conditions
%   W = zeros(size(coordinates,1),1);
%   W(unique(dirichlet)) = 0;
  
  % Solving one Newton step
  P = J\G ;
  C = C - P;
  if norm(P) < 10^(-10)
    break
  end
end


% graphic representation
figure()
subplot(1, 2, 1)
show(elements3,[],coordinates,full(C(1:M)));
title('oxygen')
subplot(1, 2, 2)
show(elements3,[],coordinates,full(C(M+1:end)));
title('carbon dioxide')