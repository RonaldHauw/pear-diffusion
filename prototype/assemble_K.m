%% Group 11 - March 12th 2020
% Assemble diffusion matrix K = [ K_u , 0 ; 0 , K_v ]

function K = assemble_K( coordinates, elements3, G2_edges, s_ur, s_vr, s_uz, s_vz, r_u, r_v )
    % coordinates       coordinates of vertices of mesh
    % elements3         index of vertices that form triangular elements
    % G2_edges          index of vertices that form the edges of outer boundary
    
    % extract useful variables
    M = size(coordinates, 1) ;
    r = coordinates(:, 1) ;
    z = coordinates(:, 2) ;

    K = zeros( 2*M, 2*M );
    for t = elements3'
                
        % area of element (can be positive or negative)
        omega = det([ ones(1,3) ; coordinates(t, :)' ]) / 2 ;
        
        % sum of r-coordinates
        sum_r = sum(coordinates(t, 1), 1) ;
        
        % for j different from i
        C_12 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(3)))*(z(t(3))-z(t(2))) ; ...
                                   (r(t(1))-r(t(3)))*(r(t(3))-r(t(2)))] ;
                               
        C_23 = 1/6 * 1/2/omega * [ (z(t(2))-z(t(1)))*(z(t(1))-z(t(3))) ; ...
                                   (r(t(2))-r(t(1)))*(r(t(1))-r(t(3)))] ;
                               
        C_13 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(2)))*(z(t(2))-z(t(3))) ; ...
                                   (r(t(1))-r(t(2)))*(r(t(2))-r(t(3)))] ;
        %
        K( t(1)  , t(2)   ) = K( t(1)  , t(2)   ) + [s_ur, s_uz] * C_12 * sum_r ;
        K( t(2)  , t(1)   ) = K( t(2)  , t(1)   ) + [s_ur, s_uz] * C_12 * sum_r ;
        %
        K( t(2)  , t(3)   ) = K( t(2)  , t(3)   ) + [s_ur, s_uz] * C_23 * sum_r ;
        K( t(3)  , t(2)   ) = K( t(3)  , t(2)   ) + [s_ur, s_uz] * C_23 * sum_r ;
        %
        K( t(1)  , t(3)   ) = K( t(1)  , t(3)   ) + [s_ur, s_uz] * C_13 * sum_r ;
        K( t(3)  , t(1)   ) = K( t(3)  , t(1)   ) + [s_ur, s_uz] * C_13 * sum_r ;

        K( t(1)+M, t(2)+M ) = K( t(1)+M, t(2)+M ) + [s_vr, s_vz] * C_12 * sum_r ;
        K( t(2)+M, t(1)+M ) = K( t(2)+M, t(1)+M ) + [s_vr, s_vz] * C_12 * sum_r ;
        %
        K( t(2)+M, t(3)+M ) = K( t(2)+M, t(3)+M ) + [s_vr, s_vz] * C_23 * sum_r ;
        K( t(3)+M, t(2)+M ) = K( t(3)+M, t(2)+M ) + [s_vr, s_vz] * C_23 * sum_r ;
        %
        K( t(1)+M, t(3)+M ) = K( t(1)+M, t(3)+M ) + [s_vr, s_vz] * C_13 * sum_r ;
        K( t(3)+M, t(1)+M ) = K( t(3)+M, t(1)+M ) + [s_vr, s_vz] * C_13 * sum_r ;
        
        % for j equal i
        C_11 = 1/6 * 1/2/omega * [ (z(t(2))-z(t(3)))^2 ; (r(t(2))-r(t(3)))^2] ;
        C_22 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(3)))^2 ; (r(t(1))-r(t(3)))^2] ;
        C_33 = 1/6 * 1/2/omega * [ (z(t(1))-z(t(2)))^2 ; (r(t(1))-r(t(2)))^2] ;
        %
        K( t(1)  , t(1)   ) = K( t(1)  , t(1)   ) + [s_ur, s_uz] * C_11 * sum_r ;
        K( t(2)  , t(2)   ) = K( t(2)  , t(2)   ) + [s_ur, s_uz] * C_22 * sum_r ;
        K( t(3)  , t(3)   ) = K( t(3)  , t(3)   ) + [s_ur, s_uz] * C_33 * sum_r ;
        %
        K( t(1)+M, t(1)+M ) = K( t(1)+M, t(1)+M ) + [s_vr, s_vz] * C_11 * sum_r ;
        K( t(2)+M, t(2)+M ) = K( t(2)+M, t(2)+M ) + [s_vr, s_vz] * C_22 * sum_r ;
        K( t(3)+M, t(3)+M ) = K( t(3)+M, t(3)+M ) + [s_vr, s_vz] * C_33 * sum_r ;

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
        K( e  , e   ) = K( e  , e   ) + r_u * [ parallel_term_1 , cross_term ; ...
                                                cross_term      , parallel_term_2 ] ;
        % in K_v
        K( e+M, e+M ) = K( e+M, e+M ) + r_v * [ parallel_term_1 , cross_term ; ...
                                                cross_term      , parallel_term_2 ] ;
    end

end