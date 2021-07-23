%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Created using Partial Differential Equation Toolbox %%%%%%%%%%%%%%%%%%%%%%%%
% Allows for the study of structural mechanics, heat transfer, and electromagnetics for an stl model%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define stl model to use
% Single stl model for complete analysis
stlLoad = 'models/ASTM.stl'; % load the stl model as a global var

% Initial idea for analyzing 'layer' by 'layer' sequentially, such that it imitates the printing process
splices = ['ASTM(1).stl','ASTM(2).stl','ASTM(3).stl', 'ASTM(4).stl'];


%% Added so user can skip GUI
prompt = 'Do you want to use the GUI or do hardcode the model data? [gui/code] \n';
guiAnswer = input(prompt,'s');

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%For  quicker testing, can skip GUI and write all characteristics here
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

strC = 'code';
if contains(guiAnswer, strC) % Alter these values if youre going to hardcode the values
    modelScale =  0.05; % Scaled down to 5% otherwise runtimes are absurd (2018 MacBook Pro)
    youngsModulus = 2100;
    poissonsRatio = 0.23333;
    massDensity = 1.2;
    constrainedFaces = 4;
    forceFace = 9;
    forceType = 'ZDisplacement';
    forceMagnitude = 1E-4;
    viewStress = 'n';
    viewStrain = 'n';
end

%Promting for user input
strG = 'gui';
if contains(guiAnswer, strG)
    prompt = {'Scale of the model (%)', 'Enter the Youngs Modulus', 'Enter the poissons ratio','Enter the mass density'};
    dlg_title = 'Define Structural Properties';
    numlines = 1;
    def = {'','','',''};
    options.Resize='on';
    options.WindowStyle='normal';
    structureProp = inputdlg(prompt, dlg_title, numlines, def, options);
end
%%
%%%%% Modeling and Meshing %%%%%%
% Mesh from stl
model = createpde(1);
importGeometry(model,stlLoad);
m = generateMesh(model, 'Hmax',5);
save('meshData.mat', 'm');
figure(6)
pdemesh(m)

%% Load/create geometry while scaling model
% Model from mesh (This is a structural model, as such, we can assign physical properties)
nodes = m.Nodes;
elements = m.Elements;
structuralmodel = createpde('structural','transient-solid'); % For stress and strain analysis
gm = importGeometry(structuralmodel, stlLoad);

if contains(guiAnswer, strC)
    scale(gm, modelScale)
else
    sc = str2double(structureProp(1));
    scale(gm, sc); % Scale down model to desired percentage
end

structuralmodel.Geometry = gm;
% structuralmodel.Geometry = geometryFromMesh(structuralmodel,nodes,elements); % (Requires very strong processing power)

%% Plot model with faces labeled
figure(7)
pdegplot(structuralmodel,'FaceLabels','on','FaceAlpha',0.5)
title('Model with faces labeled')

%%
%%%%% Evaluate structure %%%%%

%Define Structural properties based on resin/material being used
% UMA 90
if contains(guiAnswer, strC)
    structuralProperties(structuralmodel, 'YoungsModulus', youngsModulus, ... % from UMA 90 data sheet
                                          'PoissonsRatio', poissonsRatio, ... % Based on estimated transverse and axial strains 
                                          'MassDensity', massDensity); % from UMA 90 data sheet
else % Structural properties from user prompt                
    ym = str2double(structureProp(2));
    pr = str2double(structureProp(3));
    md = str2double(structureProp(4));
    structuralProperties(structuralmodel, 'YoungsModulus', ym, 'PoissonsRatio', pr, 'MassDensity', md); 
end

% Prompt user for Force constraints
if contains(guiAnswer, strG)
    prompt = {'Which faces are fixed?', 'What face(s) is the force applied to?', 'What type of force is applied?',...
        'What is the magnitude of the force?', 'Do you want to observe stress? [y/n]', 'Do you want to observe strain? [y/n]'};
    dlg_title = 'Forces';
    numlines = 1;
    def = {'','','','','',''};
    options.Resize='on';
    options.WindowStyle='normal';
    forceProp = inputdlg(prompt, dlg_title, numlines, def, options);
end

%% Define boundary constraints for constrained faces 
if contains(guiAnswer, strC)
    structuralBC(structuralmodel,'Face',constrainedFaces,'Constraint','fixed'); % Default UMA 90 ASTM example (6 left, 8 right, 4 top, 9 bottom)
else
    fc = str2double(forceProp(1));
    structuralBC(structuralmodel,'Face', fc,'Constraint','fixed'); % If testing bending forces with two or more attached faces
end

% Define Body Load (simply its body weight under gravity) Note: would be different if suspended in resin/vat
structuralBodyLoad(structuralmodel,'GravitationalAcceleration',[0;0;-9.8]); % [x,y,z]

%% Define boundary conditions for applied forces
if contains(guiAnswer, strC)
    structuralBC(structuralmodel,'Face',forceFace,forceType,forceMagnitude, 'Frequency', 50); % Nonconstant Sinusoidal displacement
    % Example Forces ------------
    % structuralBC(structuralmodel,'Face',9, 'ZDisplacement',1E-4); % Simple z-displacement
    % structuralBC(structuralmodel,'Face',9, 'YDisplacement',1E-4); % Simple y-displacement
    % structuralBC(structuralmodel,'Face',9, 'XDisplacement',1E-4); % Simple x-displacement
else

    fcc = str2double(forceProp(2));
    fct = char(forceProp(3));
    fcm =  str2double(forceProp(4));
    structuralBC(structuralmodel,'Face', fcc, fct, fcm); % Constant displacement
    % structuralBC(structuralmodel,'Face',fcc,fct,fcm, 'Frequency', 50); % Nonconstant Sinusoidal displacement
end
%% Generate mesh and 'solve' the model
generateMesh(structuralmodel,'Hmax',1); 
structuralIC(structuralmodel,'Displacement',[0;0;0],'Velocity',[0;0;0]); % Only for transient-solid
tlist = 0:0.002:0.2;
structuralResults = solve(structuralmodel,tlist);

% Evaluate reaction forces for all time steps
% reaction = evaluateReaction(structuralResults,'Face',4);

%% Stress is the force applied to a material, divided by the material’s cross-sectional area.

% Normal Stress = evaluateStress()
stressBool = "y";
if contains(guiAnswer, strG)
    if contains(forceProp(5), stressBool)
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
    end
else
    if contains(viewStress, stressBool)
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
    end
end

%% Strain is the deformation or displacement of material that results from an applied stress.

% Normal Strain - evaluateStrain()
strainBool = "y";
if contains(guiAnswer, strG)
    if contains(forceProp(6), strainBool)
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
    end
else
    if contains(viewStrain, strainBool)
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
    end
end

%% Deflection Analysis 





%$%$%$%$%$%$%$%$%$%$%$$%$%$%$%$%$%$%$%$%$%$%$%$%$$%$%$%$%$%$%$%$%$%$%
% Add prompt for deflection analysis as well as to hardcode section
%$%$%$%$%$%$%$%$%$%$%$$%$%$%$%$%$%$%$%$%$%$%$%$%$$%$%$%$%$%$%$%$%$%$%






% Import model and geometries
modelTwo = createpde('structural','static-solid');
importGeometry(modelTwo,stlLoad);

% % Define strucutural properties and boundary conditions 
structuralProperties(modelTwo,'YoungsModulus',youngsModulus, 'PoissonsRatio', poissonsRatio, 'MassDensity', massDensity); 
structuralBC(modelTwo,'Face',[14,10],'Constraint','fixed'); % Pulling down middle 'bar' for ASTM

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
% -------------------------------------------------------------------------- %

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
pdeplot3D(modelTwo,'ColorMapData',result.VonMisesStress, 'Deformation',result.Displacement, 'DeformationScaleFactor',250)
title('Simulated Model Deformation with a scale factor of %d', 350)
colormap('jet')
%}
