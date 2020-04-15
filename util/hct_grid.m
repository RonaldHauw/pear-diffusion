%% hct_grid.m provides the demanded amount of points to represent the outline
% of the border of a half-pear
%   P   :   Vector of size 1xN containing the sampled points describing the
%   boundary of the pear

function create_mesh( type, name, radius )
    
    radius = 0.02;      % 1 for example solutions, 0.02 for real tests. 
    grid_finess = 20;   % 10 for fast, 18 for accurate
    pear_shape = 0.1;   % 0.1 for pear, 1 for half circle 

    x_zero = radius-pear_shape*radius;
    a =  -x_zero * (radius^2 - x_zero^2)^(-0.5);
    y_zero = sqrt(radius^2 - x_zero^2);
    y_high = y_zero - a * x_zero;

    % % Creation of the domain
    % R1 = [ 3, 4, 0, 0, -1, -1, -1, 1, 1, -1 ]';
    % C1 = [ 1, 0, 0, radius ]';
    % P1 = [ 2, 3, 0, 0, x_zero, y_zero, y_high, y_zero ]';
    % C1 = [ C1; zeros(length(R1) - length(C1),1) ];
    % P1 = [ P1; zeros(length(R1) - length(P1),1) ];
    % % geometry description matrix
    % gd = [ R1, C1, P1 ];
    % % set formula
    % sf = 'C1-R1+P1';
    % % name space matrix that relates the colums in gd and the names in sf
    % ns = char( 'R1', 'C1', 'P1' );
    % ns = ns';

    % Creation of the domain
    R1 = [ 3, 4, 0, 0, -1, -1, -1, 1, 1, -1 ]';
    C1 = [ 1, 0, 0, radius ]';
    C2 = [ 1, 0, radius, radius/1.5 ]';
    C1 = [ C1; zeros(length(R1) - length(C1),1) ];
    C2 = [ C2; zeros(length(R1) - length(C2),1) ];
    % geometry description matrix
    gd = [ C2, C1, R1 ];
    % set formula
    sf = 'C1-R1+C2-R1';
    % name space matrix that relates the colums in gd and the names in sf
    ns = char( 'C2', 'C1', 'R1' );
    ns = ns';

    % find minimal regions that evaluate to true for the set formula sf
    dl = decsg(gd,sf,ns);

    % ivsualize mesh
    figure
    pdegplot(dl,'EdgeLabels','on','SubdomainLabels','on');
    axis equal;

    % Creation of the mesh
    model = createpde(1);
    geometryFromEdges(model,dl);
    mesh = generateMesh(model, 'GeometricOrder', 'linear', 'Hmax',radius/grid_finess*1.5,'Hmin',radius/grid_finess);
    pdeplot(model);


    %% Convert the mesh

    Nodes = [(1:size(mesh.Nodes,2))' mesh.Nodes'];
    Elements = [(1:size(mesh.Elements,2))' mesh.Elements'];
    %
    % InnerBNodes = findNodes(mesh,'region','Edge',[3 4 5]);
    % OuterBNodes = findNodes(mesh,'region','Edge',[1 6 7]);
    %
    InnerBNodes = findNodes(mesh,'region','Edge',[1 2 3]);
    OuterBNodes = findNodes(mesh,'region','Edge',[5 6 7 8]);


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

    path = '../data/meshes/' ;
    name = 'test' ;

    save( strcat(path, name), 'Nodes', 'Elements', 'InnerBEdges', 'OuterBEdges' ); 

    writematrix( Nodes,       strcat(path, name, '_Nodes.txt'),      'delimiter', 'space');
    writematrix( Elements,    strcat(path, name, '_Elements.txt'),   'delimiter', 'space');
    writematrix( InnerBEdges, strcat(path, name, '_InnerEdges.txt'), 'delimiter', 'space');
    writematrix( OuterBEdges, strcat(path, name, '_OuterEdges.txt'), 'delimiter', 'space');

end