%PEARPOINTS provides the demanded amount of points to represent the outline
% of the border of a half-pear
%   P   :   Vector of size 1xN containing the sampled points describing the
%   boundary of the pear

radius = 1.; 

x_zero = radius-0.85*radius;
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
mesh = generateMesh(model, 'GeometricOrder', 'linear', 'Hmax',radius/20+radius/20*0.5,'Hmin',radius/20);
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

save('mesh/HCTmesh3', 'Nodes', 'Elements'); 
save('mesh/HCTmesh3_Data', 'InnerBNodes', 'OuterBNodes', 'InnerBEdges', 'OuterBEdges'); 

writematrix(Nodes,'mesh/HCTmesh3_Nodes.txt','delimiter', 'space');
writematrix(Elements,'mesh/HCTmesh3_Elements.txt','delimiter', 'space');
writematrix(InnerBEdges,'mesh/HCTmesh3_InnerEdges.txt','delimiter', 'space');
writematrix(OuterBEdges,'mesh/HCTmesh3_OuterEdges.txt','delimiter', 'space');

%% Solve using C++

!cd ../; ./pear_diffusion_2 -maxit 100  -anl .5 
 
% observation: if residuals keep decreasing uniformly, the plausible
% solution is attained
% observation:

%% Plot the solution 
sol_o2 = readmatrix('mesh/solution_O_2.txt'); 
sol_co2 = readmatrix('mesh/solution_CO_2.txt');


% graphic representation
elements3   = Elements( : , 2:end ) ;
coordinates = Nodes(:, 2:end); 
figure()
subplot(1, 2, 1)
show(elements3,[],coordinates,full(sol_o2(1:end, 4)));
title('oxygen')
subplot(1, 2, 2)
show(elements3,[],coordinates,full(sol_co2(1:end, 4)));
title('carbon dioxide')



