%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Created using Partial Differential Equation Toolbox %%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Allows study of structural mechanics, heat transfer, and electromagnetics %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{%
%%%%% Modeling and Meshing %%%%%%

% load the stl model as a global var
stlLoad = 'models/ASTM(1).stl';

% Mesh from stl
model = createpde(1);
importGeometry(model,stlLoad);
m = generateMesh(model, 'Hmax',5);
save('meshData.mat', 'm');
figure(6)
pdemesh(m)

% Model from mesh (This is a structural model, as such, we can assign physical properties)
nodes = m.Nodes;
elements = m.Elements;
% structuralmodel = createpde('structural','static-solid'); % For simple static analysis 
structuralmodel = createpde('structural','transient-solid'); % For stress and strain analysis
gm = importGeometry(structuralmodel, stlLoad);
scale(gm, 0.005); % Scale down to 0.5% otherwise runtimes are absurd
structuralmodel.Geometry = gm;
% structuralmodel.Geometry = geometryFromMesh(structuralmodel,nodes,elements); % (Requires very strong processing power)

% Face label plot
figure(7)
pdegplot(structuralmodel,'FaceLabels','on','FaceAlpha',0.5)
title('Model with faces labeled')

% ---------------------------------------------------------------------------------------------------- %

%%%%% Evaluate structure %%%%%

% Define Structural properties based on resin/material being used

structuralProperties(structuralmodel, 'YoungsModulus', 2100, ... % from UMA 90 data sheet
                                      'PoissonsRatio', 0.2333333, ... % Based on estimated transverse and axial strains 
                                      'MassDensity', 1.2); % from UMA 90 data sheet
                                  
% Define boundary constraints for constrained faces (For ASTM model - 6 left, 8 right, 4 top, 9 bottom)
structuralBC(structuralmodel,'Face',4,'Constraint','fixed');
% structuralBC(structuralmodel,'Face',[6,8],'Constraint','fixed'); % If testing bending forces with two or more attached faces

% Define Body Load (simply its body weight under gravity) Note: would be different if suspended in resin/vat
structuralBodyLoad(structuralmodel,'GravitationalAcceleration',[0;0;-9.8]);

% Define boundary conditions for applied forces
% structuralBC(structuralmodel,'Face',9, 'ZDisplacement',1E-4); % Simple z-displacement
% Example Forces ------------
% structuralBC(structuralmodel,'Face',9, 'YDisplacement',1E-4); % Simple y-displacement
% structuralBC(structuralmodel,'Face',9, 'XDisplacement',1E-4); % Simple x-displacement
structuralBC(structuralmodel,'Face',9,'ZDisplacement',-1E-4, 'Frequency', 50); % Nonconstant Sinusoidal displacement

%%% Generate mesh and 'solve' the model
generateMesh(structuralmodel,'Hmax',1); 
structuralIC(structuralmodel,'Displacement',[0;0;0],'Velocity',[0;0;0]); % Only for transient-solid
tlist = 0:0.002:0.2;
structuralResults = solve(structuralmodel,tlist);

% Evaluate reaction forces for all time steps
% reaction = evaluateReaction(structuralResults,'Face',4);
%}

% ---------------------------------------------------------------------------------------------------- %

%%%% Stress is the force applied to a material, divided by the material’s cross-sectional area.

% Normal Stress = evaluateStress()
%
stress = evaluateStress(structuralResults);
figure(8)
pdeplot3D(structuralmodel,'FlowData',[stress.xx(:,end) stress.yy(:,end) stress.zz(:,end)])
title('Quiver Plot for Normal Stress at the Last Time-Step')
figure(9)
pdeplot3D(structuralmodel,'ColorMapData',stress.zz(:,end))
title('Stress in the z-direction of the Model')


% Principle Stress - evaluatePrincipalStress()
pStress = evaluatePrincipalStress(structuralResults);
I3 = pStress.s1 + pStress.s2 + pStress.s3;
I4 = pStress.s1.*pStress.s2 + pStress.s2.*pStress.s3 + pStress.s3.*pStress.s1;
% octahedral shear stress
tauOctStress = sqrt(2*(I3.^2 -3*I4))/3;
figure(10)
pdeplot3D(structuralmodel,'ColorMapData',tauOctStress(:,end))
title('Principle Stress')
%}

% ---------------------------------------------------------------------------------------------------- %

%%%% Strain is the deformation or displacement of material that results from an applied stress.

% Normal Strain - evaluateStrain()
%{
strain = evaluateStrain(structuralResults);
figure(12)
pdeplot3D(structuralmodel,'FlowData',[strain.xx(:,end) strain.yy(:,end) strain.zz(:,end)])
title('Quiver Plot for Normal Strain at the Last Time-Step')
figure(13)
pdeplot3D(structuralmodel,'ColorMapData',strain.zz(:,end))
title('Strain in z-direction of the Model')

% Priciple Strain - evaluatePrincipalStrain()
pStrain = evaluatePrincipalStrain(structuralResults);
I1 = pStrain.e1 + pStrain.e2 + pStrain.e3;
I2 = pStrain.e1.*pStrain.e2 + pStrain.e2.*pStrain.e3 + pStrain.e3.*pStrain.e1;
% octahedral shear strain
tauOctStrain = sqrt(2*(I1.^2 -3*I2))/3;
figure(14)
pdeplot3D(structuralmodel,'ColorMapData',tauOctStrain(:,end))
title('Principle Strain')
%}

% ---------------------------------------------------------------------------------------------------- %

%%%%%%% Deflection Analysis 
% % Import model and geometries
%{
modelTwo = createpde('structural','static-solid');
importGeometry(modelTwo,stlLoad);

% % Define strucutural properties and boundary conditions 
structuralProperties(modelTwo,'YoungsModulus',2100, 'PoissonsRatio', 0.2333333, 'MassDensity', 1.2); 
structuralBC(modelTwo,'Face',[14,10],'Constraint','fixed'); % Pulling down middle 'bar'

mesh = generateMesh(modelTwo); 
mv = volume(mesh)/1000;
fprintf('Volume of object = %g cm^3 \n', mv)

% % Apply forces to object
force = (mv*(1.2e-6))*9.8*0.032; % F = mgh
% structuralBoundaryLoad(modelTwo,'Face', 7,'SurfaceTraction',[0;-force;0]);
structuralBodyLoad(modelTwo, 'GravitationalAcceleration', [0;-9.8*0.005;0]);
structuralBC(modelTwo,'Face',7,'YDisplacement',-0.001); % Pulling down middle 'bar'

% ---- This would be if we wanted to apply a force at a specific vertex ---- %
% structuralBoundaryLoad(model,'Vertex',6,'Force',[0;10^4;0])

% % Solve model
result = solve(modelTwo);
minUy = min(result.Displacement.uy);
fprintf('Maximal deflection in the y-direction is %g mm. \n', minUy*100)

% % Plot desired qualities 
figure(15)
pdeplot3D(modelTwo,'ColorMapData',result.VonMisesStress)
title('Von Mises Stress')
colormap('jet')
figure(16)
pdeplot3D(modelTwo,'ColorMapData',result.Displacement.ux)
title('x-displacement')
colormap('jet')
figure(17)
pdeplot3D(modelTwo,'ColorMapData',result.Displacement.uy)
title('y-displacement')
colormap('jet')
figure(18)
pdeplot3D(modelTwo,'ColorMapData',result.Displacement.uz)
title('z-displacement')
colormap('jet')
figure(19)
pdeplot3D(modelTwo,'ColorMapData',result.VonMisesStress, 'Deformation',result.Displacement, 'DeformationScaleFactor',1)
title('Simulated Model Deformation')
colormap('jet')
%}


