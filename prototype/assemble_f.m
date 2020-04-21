%% Group 11 - April 21st 2020
%
% Assemble convection vector on boundary f = [ f_u ; f_v ]
%
% Call as
%     >> f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb )
%
% with input
%     coordinates (matrix)  r and z coordinates of each vertex
%     G2_edges    (matrix)  index of vertices that form the edges of outer boundary
%     r_u         (float)   convective mass transfer coefficient of oxygen
%     r_v         (float)   convective mass transfer coefficient of carbon dioxide
%     C_u_amb     (float)   ambient oxygen concentration
%     C_v_amb     (float)   ambient carbon dioxide concentration
%
% and output
%     f           (vector)  convection vector on boundary of finite element model

function f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb )

    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;
    
    f = zeros( 2*M, 1 ) ;
    for e = G2_edges'
        % length of edge
        len = norm(diff(coordinates(e, :)), 2) ;
        
        % terms for overlap of chapeau functions on previous and next edge
        term_1 = 1/6 * len * ( 2*r(e(1)) +   r(e(2)) ) ;
        term_2 = 1/6 * len * (   r(e(1)) + 2*r(e(2)) ) ;
        
        % in f_u
        f( e )   = f( e )   + r_u * C_u_amb * [ term_1 ; term_2 ] ;
        % in f_v
        f( M+e ) = f( M+e ) + r_v * C_v_amb * [ term_1 ; term_2 ] ;
    end
end