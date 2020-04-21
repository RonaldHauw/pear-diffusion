%% Group 11 - April 21st 2020
%
% Assemble Jacobian of nonlinearity J = dH/dC
%
% Call as
%     >> J = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
%
% with input
%     coordinates (matrix)   r and z coordinates of each vertex
%     elements3   (matrix)   index of vertices that form triangular elements
%     C           (vector)   solution concentrations at current iteration
%     dR_u_u      (function) derivative of respiration kinetics for oxygen 
%                            wrt oxygen concentration
%     dR_u_v      (function) derivative of respiration kinetics for oxygen 
%                            wrt carbon dioxide concentration
%     dR_v_u      (function) derivative of respiration kinetics for carbon 
%                            dioxide wrt oxygen concentration
%     dR_v_v      (function) derivative of respiration kinetics for carbon 
%                            dioxide wrt carbon dioxide concentration
%
% and output
%     J           (matrix)   Jacobian of nonlinearity of finite element model

function J = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v )

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;

    J = sparse( 2*M, 2*M ) ;

     % quadrature points in each vertex
%      for i = 1:M
%          % sum of areas of elements arount vertex i
%          T = elements3(any(elements3==i, 2), :) ;
%          s = 0 ;
%          for t = T'
%              s = s + abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
%          end
%          J([i, M+i],[i, M+i]) = s/3 * r(i) * [  dR_u_u(C(i), C(M+i)) ,   dR_u_v(C(i), C(M+i)) ;
%                                               - dR_v_u(C(i), C(M+i)) , - dR_v_v(C(i), C(M+i)) ] ;
%                  
%      end

   % one quadrature point in the center of element t
   for t = elements3'
        
       area = abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        
       J(t, t)     = J(t, t)     + area/9 * mean(r(t)) * dR_u_u( mean(C(t)), mean(C(M+t)) ) ;
       J(t, M+t)   = J(t, M+t)   + area/9 * mean(r(t)) * dR_u_v( mean(C(t)), mean(C(M+t)) ) ;
        
       J(M+t, t)   = J(M+t, t)   - area/9 * mean(r(t)) * dR_v_u( mean(C(t)), mean(C(M+t)) ) ;
       J(M+t, M+t) = J(M+t, M+t) - area/9 * mean(r(t)) * dR_v_v( mean(C(t)), mean(C(M+t)) ) ;
        
   end
    
end