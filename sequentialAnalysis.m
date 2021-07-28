%{
    sequentialAnalysis.m is intended to do model analysis of the supplied 
    model layer by layer to simulate forces as it prints, as opposed to the entire finished model
    like is done in main.m 


    Functions from: Sunil Bhandari (2021). slice_stl_create_path(triangles,slice_height), 
    (https://www.mathworks.com/matlabcentral/fileexchange/62113-slice_stl_create_path-triangles-slice_height), 
    MATLAB Central File Exchange. Retrieved July 23, 2021.
%}

%% Adding path to be able to call functions in folder 'upload'
addpath 'upload'
% stl_slice_and_plot() % calling example script
figure(1)
grid on
triangles = read_binary_stl_file('models/tiered v1.stl');
% triangles = read_binary_stl_file('models/turret v2.stl');
% triangles = read_binary_stl_file('models/arrowhead v1.stl');
plot_stl(triangles);
xlabel('x')
ylabel('y')
zlabel('z')

%% Read the stl file (Binary and ascii options
triangles = read_binary_stl_file('models/tiered v1.stl');
% triangles = read_binary_stl_file('models/turret v2.stl');
% triangles = read_binary_stl_file('models/arrowhead v1.stl');

% triangles = read_ascii_stl_file('models/tiered v1');

%% Rotate the triangles in desired axis,save original as well
original = triangles;
%triangles = orient_stl(triangles,'z');
triangles = rotate_stl(triangles,'x',0);

%% Designate slice increment (mm) and slice model
slice_height = 0.4;
tic;[movelist, z_slices] = slice_stl_create_path(triangles, slice_height);toc;

%% Testing something
figure(2)
hold on;
grid on
xlabel('x')
ylabel('y')
zlabel('z')
view(15,23);
axis([-Inf Inf -Inf Inf z_slices(1) z_slices(end)])
for i = 1: size(movelist,2)
    mlst_all = movelist{i};     
    if ~isempty(mlst_all)
        x = mlst_all(:,1);
        y = mlst_all(:,2);
        z = ones(size(mlst_all,1),1)*z_slices(i);
        plot3(x,y,z,'b')
    end
end   
%% Print Sliced model
figure(3)
grid on
plot_slices(movelist,z_slices, 0)
xlabel('x')
ylabel('y')
zlabel('z')

% figure(4)
% x = x(~isnan(x))';
% y = y(~isnan(y))';
% z = z(~isnan(z))';
% stlwrite('test.stl',x,y,z)



%% Old attempt at slicing stl 
%{
% Import a model and plot the mesh
importedModel = 'models/ASTM.stl';
model = createpde(1);
importGeometry(model,importedModel);
generateMesh(model, 'Hmax', 2);
% figure(1)
% pdeplot3D(model);
% view(3)

% Extract points %
vertices = model.Mesh.Nodes;
vertTransposed = transpose(vertices);
save('partialData.mat', 'vertTransposed');
x = vertTransposed(:,1);
y = vertTransposed(:,2);
z = vertTransposed(:,3);

count = 0;
for i = 1:length(x)
    if x(i) < -20
       count = count + 1;
    end
end

partial = zeros(count,3);
% Order points in ascending orer of x-axis
index = 1;
for i = 1:length(x)
    if x(i) < -20 % -20 was arbitrarily chosen to test getting partial model
        partial(index,1) = x(i);
        partial(index,2) = y(i);
        partial(index,3) = z(i);
        index = index + 1;
    end
end

% partial( ~any(partial,2), : ) = [];

%}%
xPartial = transpose(partial(:,1));
yPartial = transpose(partial(:,2));
zPartial = transpose(partial(:,3));

%%% --------------------- --------------------- --------------------- --------------------- %%%

% Plot points as a check %
figure(2)
plot3(xPartial, yPartial, zPartial, 'o', 'Color','b','MarkerSize',5,'MarkerFaceColor','#D9FFFF')

figure(3)
tri = delaunay(xPartial, yPartial);
trimesh(tri, xPartial, yPartial, zPartial);

figure(4)
shp = alphaShape(partial);
plot(shp)

figure(5)
scatter3(xPartial, yPartial,zPartial,'o')


% plot3(x,y,z,'.-')
% % Little triangles
% % The solution is to use Delaunay triangulation. Let's look at some
% % info about the "tri" variable.
% tri = delaunay(x,y);
% plot(x,y,'.')
% %
% % How many triangles are there?
% [r,c] = size(tri);
% disp(r)
% % Plot it with TRISURF
% h = trisurf(tri, x, y, z);
% axis vis3d
% % Clean it up
% axis off
% l = light('Position',[-50 -15 29]);
% set(gca,'CameraPosition',[208 -50 7687])
% lighting phong
% shading interp
% colorbar EastOutside
%}
