%% hct_grid.m provides the demanded amount of points to represent the outline
% of the border of a half-pear
%   P   :   Vector of size 1xN containing the sampled points describing the
%   boundary of the pear

function nothing = pear_grid(varargin)
    
    if nargin == 1
        grid_precision = varargin{1}; 
    else
        grid_precision = 10; 
    end

    %% Creation of the domain
    radius = .01;  % 1 for example solutions, 0.1 for real tests. 
    %grid_precision = precision; %40; % 10 for fast, 18 for accurate
    pear_height = 34.3; 
    pear_n_points = 50;

    y = linspace(0, pear_height, pear_n_points);
    x = pearpoints(y/pear_height*84.3);
    y = y./1000;
    x = x./1000;

    p = polyshape(x,y);
    t = triangulation(p);

    plot(p, 'FaceColor', 'green', 'Facealpha', 0.1, 'LineWidth', 1.5);
    hold on;
    ylabel("Heigth", 'FontSize', 30);
    xlabel("Width", 'FontSize', 30);
    axis equal
    hold off;

    %% Creation of the mesh

    model = createpde(1);
    geometryFromMesh(model,t.Points', t.ConnectivityList');
    
    mesh = generateMesh(model, 'GeometricOrder', 'linear','Hmax',radius/grid_precision*1.5,'Hmin',radius/grid_precision);
    pdeplot(model);

    clear radius grid_precision pear_height pear_n_points x y p t; 

    %% Convert the mesh

    Nodes = [(1:size(mesh.Nodes,2))' mesh.Nodes'];
    Elements = [(1:size(mesh.Elements,2))' mesh.Elements'];

    InnerBNodes = findNodes(mesh,'region','Edge',[1]);
    OuterBNodes = findNodes(mesh,'region','Edge',[2]);


    figure
    pdemesh(model,'NodeLabels','on')
    hold on
    plot(mesh.Nodes(1,InnerBNodes),mesh.Nodes(2,InnerBNodes),'or','MarkerFaceColor','g')

    figure
    pdemesh(model,'NodeLabels','on')
    hold on
    plot(mesh.Nodes(1,OuterBNodes),mesh.Nodes(2,OuterBNodes),'or','MarkerFaceColor','g')

    OuterBEdges = zeros(size(OuterBNodes, 2)-1, 3);
    for i = 1:size(OuterBNodes, 2)-1
        OuterBEdges(i, :) = [i OuterBNodes(i) OuterBNodes(i+1)]; 
    end

    InnerBEdges = zeros(size(InnerBNodes, 2)-1, 3);
    for i = 1:size(InnerBNodes, 2)-1
        InnerBEdges(i, :) = [i InnerBNodes(i) InnerBNodes(i+1)]; 
    end

    %% Export the data

    path = 'data/meshes/' ;
    name = 'pear' ;

    save( strcat(path, name), 'Nodes', 'Elements', 'InnerBEdges', 'OuterBEdges' ); 

    writematrix( Nodes,       strcat(path, name, '_Nodes.txt'),      'delimiter', 'space');
    writematrix( Elements,    strcat(path, name, '_Elements.txt'),   'delimiter', 'space');
    writematrix( InnerBEdges, strcat(path, name, '_InnerEdges.txt'), 'delimiter', 'space');
    writematrix( OuterBEdges, strcat(path, name, '_OuterEdges.txt'), 'delimiter', 'space');
    nothing = true; 
end


    %% Function describing the points along the boundary of a Conference pear

    function [y] = pearpoints(x)
    %PEARPOINTS provides the demanded amount of points to represent the outline
    % of the border of a half-pear
    %   x   :   Vector of size 1xN containing the chosen x-coordinates
    %   y   :   Vector of size 1xN containing the sampled y-coordinates
    % William Pear
    % https://www.tandfonline.com/doi/pdf/10.1080/10942912.2010.506020

    y = 4.11348 * x - 0.253106 * power(x,2) + 0.00929318 * x.^3 - 0.00019599 * x.^4 +2.08296 * 10.^(-6) * x.^5 - 8.59684 * 10.^(-9) * x.^6;

    end

    
