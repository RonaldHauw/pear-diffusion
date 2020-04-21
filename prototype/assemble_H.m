%% Group 11 - April 21st 2020
%
% Assemble nonlinearity H = [ H_u ; H_v ] with quadrature rule
%
% Call as
%     >> H = assemble_H( coordinates, elements3, C, R_u, R_v )
%
% with input
%     coordinates (matrix)   r and z coordinates of each vertex
%     elements3   (matrix)   index of vertices that form triangular elements
%     C           (vector)   solution concentrations at current iteration
%     R_u         (function) respiration kinetics for oxygen
%     R_v         (function) respiration kinetics for carbon dioxide
%
% and output
%     H           (vector)   nonlinearity of finite element model

function H = assemble_H( coordinates, elements3, C, R_u, R_v )

    % extract useful variables
    M = size(coordinates, 1) ;
        
    H = zeros( 2*M, 1 ) ;
    for t = elements3'
        % r-coordinates of triangle vertices
        r = coordinates(t, 1) ;
        
        % area of element
        area = abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        
        % quadrature points in each vertex
        % H( t )   = H( t )   + 1/3 * area .* r .* R_u(C(t), C(M+t)) ;
        % H( M+t ) = H( M+t ) - 1/3 * area .* r .* R_v(C(t), C(M+t)) ;
        
        % one quadrature point in the center of element t
        H( t )   = H( t )   + area/3 * mean(r) * R_u( mean(C(t)), mean(C(M+t)) ) ;
        H( M+t ) = H( M+t ) - area/3 * mean(r) * R_v( mean(C(t)), mean(C(M+t)) ) ;
        
    end
end