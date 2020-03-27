%% Group 11 - March 12th 2020
% build linearization of H to initialize Newton solver

function [J, l] = assemble_H_lin( coordinates, elements3, C_u_amb, C_v_amb, R_u, R_v, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form triangular elements
    
    % extract useful variables
    M = size(coordinates, 1) ;
    C_amb = [ C_u_amb*ones(M, 1) ; C_v_amb*ones(M, 1) ] ;
    
    % build linearization H(x) ~ H(C_amb) + J(C_amb)*(x-C_amb)
    H_c = assemble_H( coordinates, elements3, C_amb, R_u, R_v ) ;
    J_c = assemble_J( coordinates, elements3, C_amb, dR_u_u, dR_u_v, dR_v_u, dR_v_v ) ;
    
    J   = J_c ;
    l   = H_c - J_c*C_amb ;    
end


% function H = assemble_H_lin( coordinates, elements3, C_u_amb, C_v_amb, R_u, R_v, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
%     % coordinates       coordinates of vertices of mesh
%     % elements3         index of vertices that form triangular elements
%     
%     % extract useful variables
%     M = size(coordinates, 1) ;
%         
%     H = zeros(2*M, 1) ;
%     for t = elements3'
%         % r coordinates ot triangle vertices
%         r = coordinates(t, 1) ;
%         
%         C = 1/24 * [ 2*r(1) +   r(2) +   r(3) ; ...
%                        r(1) + 2*r(2) +   r(3) ; ...
%                        r(1) +   r(2) + 2*r(3) ] ;
%         
%         % one quadrature point in the center of element t
%         H( t )   = H( t )   + C * R_u(C_u_amb, C_v_amb) ;
%         H( M+t ) = H( M+t ) - C * R_v(C_u_amb, C_v_amb) ;
% 
%     end
% end