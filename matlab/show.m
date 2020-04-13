%% Group 11 - April 13th 2020
% Visualize two-dimensional piecewise affine function graphically
%
% Adapted from J. Alberty, C. Carstensen and S. A. Funken
% at https://www.math.hu-berlin.de/~cc/cc_homepage/
%
%    SHOW(ELEMENTS3,ELEMENTS4,COORDINATES,U) presents a two-dimensional
%    spline function graphically. ELEMENTS3 denotes a set of triangles
%    with dimension (no. of triangles) x 3 and ELEMENTS4 denotes a set of 
%    parallelograms (dimension (no. of parallelograms) x 4. Both arrays 
%    include number of nodes. The nodes have to be counted clockwise.
%    or anti-clockwise. The coordinates of the nodes are stored in an
%    (no. of coordinates) x 2 - dimensional array called COORDINATES.  
%    Its i'th row defines the x- and y- coordinate of the i'th node. U is
%    a (no. of coordinates) x 1 - dimensional array containing in the
%    i'th row the value of the spline function at the i'th node.


function show( sol, name, elements3, coordinates, c_amb )

    % extract useful parameters
    M = length(sol)/2 ;
    
    % create figure box
    figure('position', [300 100 450 400])
    
    % visualize oxygen concentration
    clim = [0, 10] ;
    
    subplot(1, 2, 1)
    box on
    hold on
    trisurf(elements3, coordinates(:,1), coordinates(:,2), full( sol(1:M) )','facecolor','interp', 'linestyle', 'none')
    hold off
    view(0,90);
    grid off
    xlim( [min(coordinates(:, 1)), max(coordinates(:, 1))] )
    ylim( [min(coordinates(:, 2)), max(coordinates(:, 2))] )
    colormap(jet(128));

    % % compute axis labels
    % labels = { 0, 0; min(clim), max(clim); c_amb(1), 'C_{amb}' } ;
    % %labels = { 0, 0; 0, 0; c_amb(1), 'C_{amb}' } ;
    % 
    % labels = sortrows(labels, 1) ;
    % 
    % if ~isempty( find( [labels{:, 2}]==c_amb(1) ) )
    %     idx = find( [labels{:, 2}]==c_amb(1) );
    %     labels(idx, :) = [] ;
    % end

    caxis(clim);
    % colorbar('YTick', [labels{:, 1}], 'YTickLabel', [labels{:, 2}], 'FontSize', 10)
    colorbar('YTick', [0:2:max(clim)], 'FontSize', 10)
    title('Oxygen [mol/m³]', 'FontSize', 10)
    
    
    % visualize carbon dioxide concentration
    clim = [0, 5] ;
    
    subplot(1, 2, 2)
    box on
    hold on
    trisurf(elements3, coordinates(:,1), coordinates(:,2), full( sol(M+1:end) )','facecolor','interp', 'linestyle', 'none')
    hold off
    view(0,90);
    grid off
    xlim( [min(coordinates(:, 1)), max(coordinates(:, 1))] )
    ylim( [min(coordinates(:, 2)), max(coordinates(:, 2))] )
    colormap(jet(128));

    % % compute axis labels
    % labels = { 0, 0; min(clim), max(clim); c_amb(1), 'C_{amb}' } ;
    % %labels = { 0, 0; 0, 0; c_amb(1), 'C_{amb}' } ;
    % 
    % labels = sortrows(labels, 1) ;
    % 
    % if ~isempty( find( [labels{:, 2}]==c_amb(1) ) )
    %     idx = find( [labels{:, 2}]==c_amb(1) );
    %     labels(idx, :) = [] ;
    % end

    caxis(clim);
    % colorbar('YTick', [labels{:, 1}], 'YTickLabel', [labels{:, 2}], 'FontSize', 10)
    colorbar('YTick', [0:2:max(clim)], 'FontSize', 10)
    title('Carbon dioxide [mol/m³]', 'FontSize', 10)
    
    sgtitle( join(['Simulated ', name]), 'FontSize', 12 )
end