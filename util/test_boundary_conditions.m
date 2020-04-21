%% Group 11 - March 24th 2020
% test correctness of boundary conditions

function test_boundary_conditions(C, G1_edges, G2_edges, coordinates, elements3, ...
                             s_ur, s_uz, s_vr, s_vz, r_u, r_v, C_u_amb, C_v_amb)

    % extract useful parameters
    % number of vertices
    M = size(coordinates, 1) ;
    % coordinates of mesh vertices
    r = coordinates(:, 1) ;
    z = coordinates(:, 2) ;
    
    
    % vector to hold residuals of boundary condition
    res_G1_o2 = zeros(size(G1_edges, 1), 1) ;
    res_G1_co2 = zeros(size(G1_edges, 1), 1) ;
    %
    res_G2_o2 = zeros(size(G2_edges, 1), 1) ;
    res_G2_co2 = zeros(size(G2_edges, 1), 1) ;

    for i = 1:size(G1_edges, 1)
        % edge to analyze
        e = G1_edges(i, :) ;
        
        % find element that has this edge
        t = elements3( any(elements3==e(1), 2) & any(elements3==e(2), 2), :) ;
        % index of third vertex in element, next to the ones that form the edge
        idx_3 = t ;
        idx_3( idx_3==e(1) | idx_3==e(2) ) = [] ;
        % solution in this element
        sol_o2  = C(t) ;
        sol_co2 = C(M+t) ;
        
        % approximate normal to boundary
        n = [0 ; -1] ;
        
        % aproximate gradient of oxygen solution
        area = abs( det([ ones(1,3) ; coordinates(t, :)' ]) ) / 2 ;
        dxi_dr   = ( z(t(3)) - z(t(1)) ) / area / 2 ;
        dxi_dz   = ( r(t(1)) - r(t(3)) ) / area / 2 ;
        deta_dr  = ( z(t(1)) - z(t(2)) ) / area / 2 ;
        deta_dz  = ( r(t(2)) - r(t(1)) ) / area / 2 ;
        %
        do2  = [ sol_o2(2) - sol_o2(1) ; sol_o2(3) - sol_o2(1) ] ;
        dco2 = [ sol_co2(2) - sol_co2(1) ; sol_co2(3) - sol_co2(1) ] ;
        %
        g_o2 = [ do2(1)*dxi_dr + do2(2)*deta_dr ; do2(1)*dxi_dz + do2(2)*deta_dz ] ;
        g_co2 = [ dco2(1)*dxi_dr + dco2(2)*deta_dr ; dco2(1)*dxi_dz + dco2(2)*deta_dz ] ;
                
        % check if boundary conditions is correct
        res_G1_o2(i)  = abs( -n' * [ s_ur, 0; 0, s_uz ] * g_o2 ) ;
        res_G1_co2(i) = abs( -n' * [ s_vr, 0; 0, s_vz ] * g_co2 ) ;
        
    end

    for i = 1:size(G1_edges, 1)
        % edge to analyze
        e = G1_edges(i, :) ;
        
        % find element that has this edge
        t = elements3( any(elements3==e(1), 2) & any(elements3==e(2), 2), :) ;
        % index of third vertex in element, next to the ones that form the edge
        idx_3 = t ;
        idx_3( idx_3==e(1) | idx_3==e(2) ) = [] ;
        % solution in this element
        sol_o2  = C(t) ;
        sol_co2 = C(M+t) ;
        
        % approximate normal to boundary
        dr = diff(r(e)) ;
        dz = diff(z(e)) ;
        m  = [mean(r(e)) ; mean(z(e))] ;
        % compute the two vertors normal to this edge
        n1 = [dr ; -dz]/norm([dr ; -dz]) ;
        n2 = [-dr ; dz]/norm([-dr ; dz]) ;
        % select the one that points outward the domain 
        if norm( (n1+m) - coordinates(idx_3, :)' ) < norm( (n2+m) - coordinates(idx_3, :)' )
            n = n2 / norm(n2) ;
        else
            n = n1 / norm(n1) ;
        end
        
        % aproximate gradient of oxygen solution
        area = abs( det([ ones(1,3) ; coordinates(t, :)' ]) ) / 2 ;
        dxi_dr   = ( z(t(3)) - z(t(1)) ) / area / 2 ;
        dxi_dz   = ( r(t(1)) - r(t(3)) ) / area / 2 ;
        deta_dr  = ( z(t(1)) - z(t(2)) ) / area / 2 ;
        deta_dz  = ( r(t(2)) - r(t(1)) ) / area / 2 ;
        %
        do2  = [ sol_o2(2) - sol_o2(1) ; sol_o2(3) - sol_o2(1) ] ;
        dco2 = [ sol_co2(2) - sol_co2(1) ; sol_co2(3) - sol_co2(1) ] ;
        %
        g_o2 = [ do2(1)*dxi_dr + do2(2)*deta_dr ; do2(1)*dxi_dz + do2(2)*deta_dz ] ;
        g_co2 = [ dco2(1)*dxi_dr + dco2(2)*deta_dr ; dco2(1)*dxi_dz + dco2(2)*deta_dz ] ;
                
        % check if boundary conditions is correct
        res_G2_o2(i) = abs( -n' * [ s_ur, 0; 0, s_uz ] * g_o2 - r_u*(mean(sol_o2)-C_u_amb) ) ;
        res_G2_co2(i) = abs( -n' * [ s_vr, 0; 0, s_vz ] * g_co2 - r_v*(mean(sol_co2)-C_v_amb) ) ;
                
    end

    
    figure
    % plot residuals on inner boundary    
    subplot(2, 2, 1)
    semilogy(res_G1_o2)
    xlim([1, length(res_G1_o2)])
    xlabel('Edge label') ;
    title('Oxygen inner') ;
    
    subplot(2, 2, 2)
    semilogy(res_G1_co2)
    xlim([1, length(res_G1_co2)])
    xlabel('Edge label') ;
    title('Carbon dioxide inner') ;
    
    sgtitle('Residual on boundary condition') ;
    
    
    % plot residuals on outer boundary
    
    subplot(2, 2, 3)
    semilogy(res_G2_o2)
    xlim([1, length(res_G2_o2)])
    xlabel('Edge label') ;
    title('Oxygen outer') ;
    
    subplot(2, 2, 4)
    semilogy(res_G2_co2)
    xlim([1, length(res_G2_co2)])
    xlabel('Edge label') ;
    title('Carbon dioxide outer') ;
    
    sgtitle('Residual on boundary condition') ;
end