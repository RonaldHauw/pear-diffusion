%% Group 11 - April 15th 2020
%
% Generate pear mesh (nodes, elements, inner boundary and outer boundary).
%
% Call as
%     >> create_mesh( type, file_name, radius, granularity )
% with
%     type        (string)   shape of mesh to generate
%                            for a circle type = 'circle', 'c' or 'hc'
%                            for a circle with triangle type = 'circle triangle', 'ct' or 'hct'
%                            for a two circles type = 'double circle', 'dc' or 'cc'
%     name        (string)   name of file to save mesh (in folder data/meshes/)
%     radius      (float)    radius of pear in meters
%     granularity (int)      how fine should the regularization be
%                            5 is rough, 15 is fine, 30 is very fine

function create_mesh( type, name, radius, granularity )
    
    switch lower(type)
        
        case {'circle', 'c', 'hc'}
            
            % Creation of the domain
            R1 = [ 3, 4, 0, 0, -1, -1, -1, 1, 1, -1 ]';
            C1 = [ 1, 0, 0, radius ]';
            C1 = [ C1; zeros(length(R1) - length(C1),1) ];
            % geometry description matrix
            gd = [ R1, C1 ];
            % set formula
            sf = 'C1-R1';
            % name space matrix that relates the colums in gd and the names in sf
            ns = char( 'R1', 'C1' );
            
            in_edge  = [1] ;
            out_edge = [2 3] ;
            
        case {'circle triangle', 'ct', 'hct'}
            
            pear_shape = 0.1;
            x_zero = radius-pear_shape*radius;
            a      =  -x_zero * (radius^2 - x_zero^2)^(-0.5);
            y_zero = sqrt(radius^2 - x_zero^2);
            y_high = y_zero - a * x_zero;

            % Creation of the domain
            R1 = [ 3, 4, 0, 0, -1, -1, -1, 1, 1, -1 ]';
            C1 = [ 1, 0, 0, radius ]';
            P1 = [ 2, 3, 0, 0, x_zero, y_zero, y_high, y_zero ]';
            C1 = [ C1; zeros(length(R1) - length(C1),1) ];
            P1 = [ P1; zeros(length(R1) - length(P1),1) ];
            % geometry description matrix
            gd = [ R1, C1, P1 ];
            % set formula
            sf = 'C1-R1+P1';
            % name space matrix that relates the colums in gd and the names in sf
            ns = char( 'R1', 'C1', 'P1' );
            
            in_edge  = [3 4 5] ;
            out_edge = [1 6 7] ;
        
        case {'double circle', 'dc', 'cc'}
            
            % Creation of the domain
            R1 = [ 3, 4, 0, 0, -1, -1, -1, 1, 1, -1 ]';
            C1 = [ 1, 0, 0, radius ]';
            C2 = [ 1, 0, 1.2*radius, radius/1.5 ]';
            C1 = [ C1; zeros(length(R1) - length(C1),1) ];
            C2 = [ C2; zeros(length(R1) - length(C2),1) ];
            % geometry description matrix
            gd = [ C2, C1, R1 ];
            % set formula
            sf = 'C1-R1+C2-R1';
            % name space matrix that relates the colums in gd and the names in sf
            ns = char( 'C2', 'C1', 'R1' );
            
            in_edge  = [1 2 3] ;
            out_edge = [5 6 7 8] ;
        
        otherwise
            error("Did not understand which mesh to generate. Run help help create_mesh for more information")
    end
    
    ns = ns';
    % find minimal regions that evaluate to true for the set formula sf
    dl = decsg(gd,sf,ns);

    % visualize mesh
    figure('position', [100, 100, 800, 300])
    subplot(1, 3, 1)
    pdegplot(dl,'EdgeLabels','on','SubdomainLabels','on');
    axis equal;
    title("Mesh")

    % Creation of the mesh
    model = createpde(1);
    geometryFromEdges(model,dl);
    mesh = generateMesh(model, 'GeometricOrder', 'linear', 'Hmax',radius/granularity*1.5,'Hmin',radius/granularity);
    pdeplot(model);


    %% Convert the mesh

    Nodes = [(1:size(mesh.Nodes,2))' mesh.Nodes'];
    Elements = [(1:size(mesh.Elements,2))' mesh.Elements'];
    %
    InnerBNodes = findNodes(mesh,'region','Edge', in_edge );
    OuterBNodes = findNodes(mesh,'region','Edge', out_edge);

    subplot(1, 3, 2)
    pdemesh(model,'NodeLabels','on')
    hold on
    plot(mesh.Nodes(1,InnerBNodes),mesh.Nodes(2,InnerBNodes),'or','MarkerFaceColor','g')
    title("Inner boundary nodes")
    
    subplot(1, 3, 3)
    pdemesh(model,'NodeLabels','on')
    hold on
    plot(mesh.Nodes(1,OuterBNodes),mesh.Nodes(2,OuterBNodes),'or','MarkerFaceColor','g')
    title("Outer boundary nodes")
    
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

    save( strcat(path, name), 'Nodes', 'Elements', 'InnerBEdges', 'OuterBEdges' ); 

    writematrix( Nodes,       strcat(path, name, '_Nodes.txt'),      'delimiter', 'space');
    writematrix( Elements,    strcat(path, name, '_Elements.txt'),   'delimiter', 'space');
    writematrix( InnerBEdges, strcat(path, name, '_InnerEdges.txt'), 'delimiter', 'space');
    writematrix( OuterBEdges, strcat(path, name, '_OuterEdges.txt'), 'delimiter', 'space');

end