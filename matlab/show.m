function show(elements3,elements4,coordinates,u, clim, c_amb)
%SHOW   Presents two-dimensional piecewise affine function graphically.
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

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <show.m> in $(HOME)/acf/fem2d/ and
%                  in $(HOME)/acf/fem2d_heat/ and
%                  in $(HOME)/acf/fem2d_nonlinear/.

box on
hold on
trisurf(elements3,coordinates(:,1),coordinates(:,2),u','facecolor','interp', 'linestyle', 'none')
trisurf(elements4,coordinates(:,1),coordinates(:,2),u','facecolor','interp', 'linestyle', 'none')
hold off
view(0,90);
% axis equal
grid off
xlim( [min(coordinates(:, 1)), max(coordinates(:, 1))] )
ylim( [min(coordinates(:, 2)), max(coordinates(:, 2))] )
title('Solution of the Problem')
colormap(jet(128));

% % compute axis labels
% labels = { 0, 0; min(clim), max(clim); c_amb, 'C_{amb}' } ;
% %labels = { 0, 0; 0, 0; c_amb, 'C_{amb}' } ;
% 
% labels = sortrows(labels, 1) ;
% 
% if ~isempty( find( [labels{:, 2}]==c_amb ) )
%     idx = find( [labels{:, 2}]==c_amb );
%     labels(idx, :) = [] ;
% end

caxis(clim);
% colorbar('YTick', [labels{:, 1}], 'YTickLabel', [labels{:, 2}], 'FontSize', 10)
colorbar('YTick', [0:2:max(clim)], 'FontSize', 10)

end