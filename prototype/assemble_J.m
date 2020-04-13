%% Group 11 - March 27th 2020
% Assemble Jacobian J = dH/dC

function J = assemble_J( coordinates, elements3, C, dR_u_u, dR_u_v, dR_v_u, dR_v_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form trinagular elements

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;

    J = zeros( 2*M, 2*M ) ;

%     % quadrature points in each vertex
%     for i = 1:M
%         % sum of areas of elements arount vertex i
%         T = elements3(any(elements3==i, 2), :) ;
%         s = 0 ;
%         for t = T'
%             s = s + abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
%         end
%         J([i, M+i],[i, M+i]) = s/3 * r(i) * [  dR_u_u(C(i), C(M+i)) ,   dR_u_v(C(i), C(M+i)) ;
%                                              - dR_v_u(C(i), C(M+i)) , - dR_v_v(C(i), C(M+i)) ] ;
%                 
%     end

    % one quadrature point in the center of element t
    for t = elements3'
        
        area = abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        
        J(t, t)     = J(t, t)     + area/9 * mean(r(t)) * dR_u_u( mean(C(t)), mean(C(M+t)) ) ;
        J(t, M+t)   = J(t, M+t)   + area/9 * mean(r(t)) * dR_u_v( mean(C(t)), mean(C(M+t)) ) ;
        
        J(M+t, t)   = J(M+t, t)   - area/9 * mean(r(t)) * dR_v_u( mean(C(t)), mean(C(M+t)) ) ;
        J(M+t, M+t) = J(M+t, M+t) - area/9 * mean(r(t)) * dR_v_v( mean(C(t)), mean(C(M+t)) ) ;
        
    end
end