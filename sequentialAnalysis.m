%{
sequentialAnalysis.m is intended to do model analysis of the supplied 
model layer by layer as it prints, as opposed to the entire finished model
like is done pProc.m 
%}

% Import a model and plot the mesh
importedModel = 'ASTM.stl';
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
