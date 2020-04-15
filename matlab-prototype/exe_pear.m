%EXE_PEAR 
%   Initialises the mesh, exports it to .txt format, runs the C++ program
%   and then visualises the output

%   Approximates the outline of a half-pear with a polynomial of the sixth
%   order


%% Creation of the domain

radius = .01;  % 1 for example solutions, 0.1 for real tests. 
grid_precision = 3; % 10 for fast, 18 for accurate
pear_height = 84.3; 
pear_n_points = 30;

y = linspace(0, pear_height, pear_n_points);
x = pearpoints(y);

y = y./2000;
x = x./2000;

% [x,y] = click;
% x = x/7; y = y/7;
% x = x - x(1); y(end) = y(end-1); y = y-y(end);

p = polyshape(x,y);
t = triangulation(p);

plot(p, 'FaceColor', 'green', 'Facealpha', 0.1, 'LineWidth', 1.5);
hold on;
ylabel("Heigth (mm)", 'FontSize', 30);
xlabel("Width (mm)", 'FontSize', 30);
% xlim([0 40]);
% ylim([0 90]);
axis equal
hold off;

%% Creation of the mesh

model = createpde(1);
geometryFromMesh(model,t.Points', t.ConnectivityList');
mesh = generateMesh(model, 'GeometricOrder', 'linear','Hmax',0.0015,'Hmin',0.0008);
pdeplot(model);
    
clear radius grid_precision pear_height pear_n_points x y p t; 

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

clear InnerBEdges InnerBNodes OuterBEdges OuterBNodes; 

%% Solve using C++

!cd ../; ./pear_diffusion_2 -maxit 100  -anl 0.05
 

%% Plot the solution 
sol_o2 = readmatrix('mesh/solution_O_2.txt'); 
sol_co2 = readmatrix('mesh/solution_CO_2.txt');


% Graphic representation
elements3   = Elements( : , 2:end ) ;
coordinates = Nodes(:, 2:end); 
figure()
subplot(1, 2, 1)
show(elements3,[],coordinates,full(sol_o2(1:end, 4)));
title('oxygen')
subplot(1, 2, 2)
show(elements3,[],coordinates,full(sol_co2(1:end, 4)));
title('carbon dioxide')

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

function [x,y] = click

	% initialiseer figuur

    img = imread('pear.png');
    
	figure(1); clf
    hold on;
    imagesc([0 0.5], [0 1], flipud(img));
	axis([0 0.5 0 1]);
	axis equal
	title('Click left to draw polyline, click right to terminate')
	hold on;

	% herhaal tot andere dan linkermuisknop ingedrukt
	x = []; y = [];
	while(1)
		[px,py,button] = ginput(1);
		if( button ~= 1 )
			break;
		else
			x = [x px] ; y = [y py];
			if( length(x) > 1 )
				plot(x([end-1 end]),y([end-1 end]),'b-');
			end
			plot(px,py,'+');
		end
	end

	hold off;
end



