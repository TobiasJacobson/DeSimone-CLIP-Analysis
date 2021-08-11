%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Created using Partial Differential Equation Toolbox %%%%%%%%%%%%%%%%%%%%%%%%
% Allows for the study of structural mechanics, heat transfer, and electromagnetics for an stl model%
%%% Note:  In order to run the script, you must install the Partial Differential Equation Toolbox %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Issue: Won't apply force to all defined needle points %%%




%{
function [] = singleModelAnalysis(stl, ms, ym, pr, md, cf, ff, ft, fm, vstress, vstrain, dd, ds, daf, dff, dft, dfm)
%% Define stl model to use
% Single stl model for complete analysis
stlLoad = stl; % load the stl model as a global var


%%%%%%%%%%
% Used for stress and strain
% At the moment these values are for the CARBON UMA 90 resin
%%%%%%%%%%
modelScale = ms; % Model scale: Scaled down to 5% otherwise runtimes are in excess of 60min (2018 MacBook Pro)
youngsModulus = ym; % From uma 90 data sheet
poissonsRatio = pr; % Based on estimated transverse and axial strains 
massDensity = md; % From uma 90 data sheet
constrainedFaces = cf; % Face(s) attached to the baseplate or that are simply fixed in place
forceFace = ff; % Face(s) which have the force applied to them
forceType = ft'; % Options: Displacement [x;y;z], XDisplacement, YDisplacement, ZDisplacement, RDisplacement
forceMagnitude = fm;  % The displacement force magnitude
viewStress = vstress; % View stress and principle stress analysis? [y/n]
viewStrain = vstrain; % View strain and principle strain analysis? [y/n]

% Used for deflection analysis
defHBool = dd; % Do deflection analysis [y/n]
defScale = ds; % Model scale: Scaled down to 5% otherwise runtimes are in excess of 60min (2018 MacBook Pro)
defAFace = daf; % Face(s) attached to the baseplate or that are simply fixed in place
defFFace = dff; % Face(s) which have the force applied to them
defFType = dft; % Options: Displacement [x;y;z], XDisplacement, YDisplacement, ZDisplacement, RDisplacement
defFMag = dfm; % The displacement force magnitude

%}


%%%%% -------- %%%%%%% This is just being used for debugging purposes %%%%%%% -------- %%%%%

stlLoad = 'models/tiered v1.stl';
modelScale =  0.05; % Model scale: Scaled down to 5% otherwise runtimes are in excess of 60min (2018 MacBook Pro)
youngsModulus = 2100; % From uma 90 data sheet
poissonsRatio = 0.23333; % Based on estimated transverse and axial strains 
massDensity = 1.2; % From uma 90 data sheet
constrainedFaces = [1, 10, 11, 20, 21, 30, 36, 41, 46, 51, 56, 61, 66, 67, 72, 77, 82, 87, 92, ... 
    101, 106, 111, 116, 121, 126, 127, 132, 137, 142, 147, 152, 161, 166, 171, 176, ...
    181, 186, 187, 192, 197, 202, 207, 212]; % Face(s) attached to the baseplate or that are simply fixed in place
forceFace = 9; % Face(s) which have the force applied to them
forceVertex = [5, 6, 15, 16, 25, 26, 39, 44, 49, 54, 59, 64, 73, 78, 83, 88, 93, 98, ...
        99, 104, 109, 114, 119, 124, 133, 138, 143, 148, 153, 158, 159, 164, 169, ...
        174, 179, 184, 193, 198, 203, 208, 213, 218];
forceType = 'ZDisplacement'; % Options: Displacement [x;y;z], XDisplacement, YDisplacement, ZDisplacement, RDisplacement
forceMagnitude = -10;  % The displacement force magnitude
viewStress = 'n'; % View stress and principle stress analysis? [y/n]
viewStrain = 'n'; % View strain and principle strain analysis? [y/n]

% Used for deflection analysis
defHBool = 'y'; % Do deflection analysis [y/n]
defScale = 0.05; % Model scale: Scaled down to 5% otherwise runtimes are in excess of 60min (2018 MacBook Pro)
defAFace = [1, 10, 11, 20, 21, 30, 36, 41, 46, 51, 56, 61, 66, 67, 72, 77, 82, 87, 92, ... 
    101, 106, 111, 116, 121, 126, 127, 132, 137, 142, 147, 152, 161, 166, 171, 176, ...
    181, 186, 187, 192, 197, 202, 207, 212]; % Face(s) attached to the baseplate or that are simply fixed in place
defFFace = 7; % Face(s) which have the force applied to them
defFVertex = [5, 6, 15, 16, 25, 26, 39, 44, 49, 54, 59, 64, 73, 78, 83, 88, 93, 98, ...
        99, 104, 109, 114, 119, 124, 133, 138, 143, 148, 153, 158, 159, 164, 169, ...
        174, 179, 184, 193, 198, 203, 208, 213, 218];
defFType = 'ZDisplacement'; % Options: Displacement [x;y;z], XDisplacement, YDisplacement, ZDisplacement, RDisplacement
defFMag = -0.0001; % The displacement force magnitude

%%%%%%% -------- %%%%%%% -------- %%%%%%% -------- %%%%%%% -------- %%%%%%% -------- %%%%%%%
%%
%%%%% PDE Meshing %%%%%%
model = createpde();
importGeometry(model, stlLoad);
m = generateMesh(model, 'Hmax',0.5); % Generate mesh from stl geometry
save('meshData.mat', 'm'); % save mesh data incase needed for other use 
figure(1)
pdemesh(m)

%% Load/create geometry while scaling model
% Model from mesh (This is a structural model, as such, we can assign physical properties)
nodes = m.Nodes;
elements = m.Elements;
structuralmodel = createpde('structural','static-solid'); % Define PDE model type
gm = importGeometry(structuralmodel, stlLoad);
scale(gm, modelScale)
structuralmodel.Geometry = gm;
% structuralmodel.Geometry = geometryFromMesh(structuralmodel,nodes,elements); % (Requires very strong processing power)

%% Plot model with faces/vertices labeled
figure(2)
%
% For vertices 
pdegplot(structuralmodel,'VertexLabels','on','FaceAlpha',0.5) % 50% transparency 
title('Model with vertices labeled')
%}

%{
% For faces
pdegplot(structuralmodel,'FaceLabels','on','FaceAlpha',0.5) % 50% transparency 
title('Model with faces labeled')
%}

%%
%%%%% Evaluate structure based on forces applied %%%%%

%Define Structural properties based on resin/material being used
structuralProperties(structuralmodel, 'YoungsModulus', youngsModulus, ... 
                                      'PoissonsRatio', poissonsRatio, ... 
                                      'MassDensity', massDensity); 

%% Define boundary constraints for constrained faces 
structuralBC(structuralmodel,'Face',constrainedFaces,'Constraint','fixed'); %  UMA 90 ASTM example (6 left, 8 right, 4 top, 9 bottom)
% Define Body Load (simply its body weight under gravity) Note: would be different if suspended in resin/vat
structuralBodyLoad(structuralmodel,'GravitationalAcceleration',[0;0;-9.8]); % [x,y,z]

%% Define boundary conditions for applied forces
% For faces -- 
% structuralBC(structuralmodel,'Face',forceFace,forceType,forceMagnitude, 'Frequency', 50); % Nonconstant Sinusoidal displacement

% For vertices -- 
structuralBC(structuralmodel,'Vertex',forceVertex,forceType,forceMagnitude, 'Frequency', 50);

% Other Example Forces ------------
% structuralBC(structuralmodel,'Face',9, 'ZDisplacement',1E-4); % Constant z-displacement
% structuralBC(structuralmodel,'Face',9, 'YDisplacement',1E-4); % Constant y-displacement
% structuralBC(structuralmodel,'Face',9, 'XDisplacement',1E-4); % Constant x-displacement

%% Generate mesh and 'solve' the model
generateMesh(structuralmodel,'Hmax',0.5); 
tlist = 0:0.002:0.2;
structuralResults = solve(structuralmodel,tlist);

%% Stress is the force applied to a material divided by the material’s cross-sectional area.

% Normal Stress = evaluateStress()
stressBool = "y";
if contains(viewStress, stressBool)
    % Stress - evaluateStress()
    stress = evaluateStress(structuralResults);
    figure(3)
    pdeplot3D(structuralmodel,'FlowData',[stress.xx(:,end) stress.yy(:,end) stress.zz(:,end)])
    title('Quiver Plot for Normal Stress at the Last Time-Step')
    figure(4)
    pdeplot3D(structuralmodel,'ColorMapData',stress.zz(:,end))
    title('Stress in the z-direction of the Model')

    % Principle Stress - evaluatePrincipalStress()
    pStress = evaluatePrincipalStress(structuralResults);
    I3 = pStress.s1 + pStress.s2 + pStress.s3;
    I4 = pStress.s1.*pStress.s2 + pStress.s2.*pStress.s3 + pStress.s3.*pStress.s1;
    % octahedral shear stress
    tauOctStress = sqrt(2*(I3.^2 -3*I4))/3;
    figure(5)
    pdeplot3D(structuralmodel,'ColorMapData',tauOctStress(:,end))
    title('Principle Stress')
end

%% Strain is the deformation or displacement of material that results from an applied stress.

% Normal Strain - evaluateStrain()
strainBool = "y";
if contains(viewStrain, strainBool)
    % Strain - evaluateStrain()
    strain = evaluateStrain(structuralResults);
    figure(6)
    pdeplot3D(structuralmodel,'FlowData',[strain.xx(:,end) strain.yy(:,end) strain.zz(:,end)])
    title('Quiver Plot for Normal Strain at the Last Time-Step')
    figure(7)
    pdeplot3D(structuralmodel,'ColorMapData',strain.zz(:,end))
    title('Strain in z-direction of the Model')

    % Priciple Strain - evaluatePrincipalStrain()
    pStrain = evaluatePrincipalStrain(structuralResults);
    I1 = pStrain.e1 + pStrain.e2 + pStrain.e3;
    I2 = pStrain.e1.*pStrain.e2 + pStrain.e2.*pStrain.e3 + pStrain.e3.*pStrain.e1;
    % octahedral shear strain
    tauOctStrain = sqrt(2*(I1.^2 -3*I2))/3;
    figure(8)
    pdeplot3D(structuralmodel,'ColorMapData',tauOctStrain(:,end))
    title('Principle Strain')
end

%% Deflection Analysis

isYes = 'y';
if contains(defHBool, isYes)
    % Import model and geometries
    modelTwo = createpde('structural','static-solid');
    importGeometry(modelTwo,stlLoad);
    
    figure(99)
    pdegplot(modelTwo,'FaceLabels','on','FaceAlpha',0.5) % 50% transparency 
    title('Model with faces labeled')

    structuralProperties(modelTwo,'YoungsModulus',youngsModulus, 'PoissonsRatio', poissonsRatio, 'MassDensity', massDensity); 
    structuralBC(modelTwo,'Face',defAFace,'Constraint','fixed');

    mesh = generateMesh(modelTwo); 
    mv = volume(mesh)/1000;
    fprintf('Volume of object = %g cm^3 \n', mv)

    % Apply forces to object
	
%     force = -9.8*defScale;
%     structuralBodyLoad(modelTwo, 'GravitationalAcceleration', [0;0;force]); % Force of gravity acting on object

    % ----- This would be if we wanted to apply a force at a specific face(s) ----- %
%     structuralBC(modelTwo,'Face',defFFace,defFType,defFMag); 
    % ----------------------------------------------------------------------------- %
    
    % ---- This would be if we wanted to apply a force at a specific vertex(s) ---- %
    structuralBoundaryLoad(modelTwo,'Vertex', defFVertex ,'Force',[0;0;1])
    % ----------------------------------------------------------------------------- %

    % Solve model
    result = solve(modelTwo);
    minUx = min(result.Displacement.ux);
    minUy = min(result.Displacement.uy);
    minUz = min(result.Displacement.uz);
    fprintf('Maximal deflection in the x-direction is %g mm. \n', minUx*100)
    fprintf('Maximal deflection in the y-direction is %g mm. \n', minUy*100)
    fprintf('Maximal deflection in the z-direction is %g mm. \n', minUz*100)

    % Plot desired qualities 
    figure(9)
    pdeplot3D(modelTwo,'ColorMapData',result.VonMisesStress)
    title('Von Mises Stress')
    colormap('jet')
    figure(10)
    pdeplot3D(modelTwo,'ColorMapData',result.Displacement.ux)
    title('x-displacement')
    colormap('jet')
    figure(11)
    pdeplot3D(modelTwo,'ColorMapData',result.Displacement.uy)
    title('y-displacement')
    colormap('jet')
    figure(12)
    pdeplot3D(modelTwo,'ColorMapData',result.Displacement.uz)
    title('z-displacement')
    colormap('jet')
    figure(13)
    pdeplot3D(modelTwo,'ColorMapData',result.VonMisesStress, 'Deformation',result.Displacement, 'DeformationScaleFactor',0)
    title('Simulated Model Deformation with a scale factor of 0%')
    colormap('jet')
end
% end
%% End Code
