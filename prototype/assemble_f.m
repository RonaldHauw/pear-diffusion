%% Group 11 - March 12th 2020
% Assemble f = [ f_u ; f_v ]

function f = assemble_f( coordinates, G2_edges, r_u, r_v, C_u_amb, C_v_amb )
    % coordinates       coordinates of vertices of mesh
    % G2_edges          index of vertices that form the edges of outer boundary 

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
        f( e )   = f( e )   + r_u * C_u_amb * [ term_1 ; term_2 ] ;
        % in f_v
        f( M+e ) = f( M+e ) + r_v * C_v_amb * [ term_1 ; term_2 ] ;
    end
end