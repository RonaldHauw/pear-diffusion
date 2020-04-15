%PEARPOINTS provides the demanded amount of points to represent the outline
% of the border of a half-pear
%   P   :   Vector of size 1xN containing the sampled points describing the
%   boundary of the pear

radius = 0.02;  % 1 for example solutions, 0.1 for real tests. 
grid_finess = 30; % 10 for fast, 18 for accurate
pear_shape = 0.1; % 0.1 for pear, % 0.85 for half circle 

x_zero = radius-pear_shape*radius;
a = - x_zero * (radius^2 - x_zero^2)^(-0.5);
y_zero = sqrt(radius^2-x_zero^2);
y_high = y_zero - a * x_zero;

% Creation of the domain
R1 = [3,4,0, 0, -1, -1, -1 ,1, 1, -1]';
C1 = [1,0,0,radius]';
P1 = [2, 3, 0, 0, x_zero, y_zero, y_high, y_zero]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
P1 = [P1;zeros(length(R1) - length(P1),1)];
gd = [R1,C1, P1];
sf = 'C1-R1+P1';
ns = char('R1','C1','P1');
ns = ns';
dl = decsg(gd,sf,ns);
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
InnerBNodes = findNodes(mesh,'region','Edge',[3 4 5]);
OuterBNodes = findNodes(mesh,'region','Edge',[1 6 7]);

%figure
%pdemesh(model,'NodeLabels','on')
%hold on
%plot(mesh.Nodes(1,InnerBNodes),mesh.Nodes(2,InnerBNodes),'or','MarkerFaceColor','g')

%figure
%pdemesh(model,'NodeLabels','on')
%hold on
%plot(mesh.Nodes(1,OuterBNodes),mesh.Nodes(2,OuterBNodes),'or','MarkerFaceColor','g')

OuterBEdges = zeros(size(OuterBNodes, 2)-1, 3);
for i = 1:size(OuterBNodes, 2)-1
    OuterBEdges(i, :) = [i OuterBNodes(i) OuterBNodes(i+1)]; 
end

InnerBEdges = zeros(size(InnerBNodes, 2)-1, 3);
for i = 1:size(InnerBNodes, 2)-1
    InnerBEdges(i, :) = [i InnerBNodes(i) InnerBNodes(i+1)]; 
end

%% Export the data

save('../data/meshes/HCTmesh3', 'Nodes', 'Elements', 'InnerBEdges', 'OuterBEdges'); 

writematrix(Nodes,'../data/meshes/HCTmesh3_Nodes.txt','delimiter', 'space');
writematrix(Elements,'../data/meshes/HCTmesh3_Elements.txt','delimiter', 'space');
writematrix(InnerBEdges,'../data/meshes/HCTmesh3_InnerEdges.txt','delimiter', 'space');
writematrix(OuterBEdges,'../data/meshes/HCTmesh3_OuterEdges.txt','delimiter', 'space');