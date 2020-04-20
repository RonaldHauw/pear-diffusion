%% Group 11 - March 27th 2020
% Assemble H = [ H_u ; H_v]

function H = assemble_H( coordinates, elements3, C, R_u, R_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form triangular elements

    % extract useful variables
    M = size(coordinates, 1) ;
        
    H = sparse( 2*M, 1 ) ;
    for t = elements3'
        % r-coordinates of triangle vertices
        r = coordinates(t, 1) ;
        
        % area of element
        area = abs(det([ ones(1,3) ; coordinates(t, :)' ])) / 2 ;
        
%         % quadrature points in each vertex
         H( t )   = H( t )   + 1/3 * area .* r .* R_u(C(t), C(M+t)) ;
         H( M+t ) = H( M+t ) - 1/3 * area .* r .* R_v(C(t), C(M+t)) ;
        
        % one quadrature point in the center of element t
%        H( t )   = H( t )   + area/3 * mean(r) * R_u( mean(C(t)), mean(C(M+t)) ) ;
%        H( M+t ) = H( M+t ) - area/3 * mean(r) * R_v( mean(C(t)), mean(C(M+t)) ) ;
        
    end
end