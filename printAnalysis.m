%% 
%{
    printAnalysis.m is intented to run singleModelAnalysis at each layer
    of a print. Showing model/print characteristics as the print progresses

    Must first define model properties and material properties
    The current proerties are for an the tiered v1 micro needle and the UMA 90 @CARBON resin

%}

%%
% Define resin properties and model characteristics
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
viewDeflection = 'y'; % Do deflection analysis [y/n]
defScale = 0; %0.01 % Deflection Deformation Scale


singleModelAnalysis(stlLoad, modelScale, youngsModulus, poissonsRatio, massDensity, constrainedFaces, forceFace, forceVertex, ...
    forceType, forceMagnitude, viewStress, viewStrain, viewDeflection, defScale)

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Call something from sequential analysis to return stl's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call here 
stlModel = output from sequential analysis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run all stls through singleModelAnalysis 
numStls = 4;
for i = 1:numStls
    
    singleModelAnalysis(stlModel, modelScale, youngsModulus, poissonsRatio, massDensity, constrainedFaces, forceFace, forceType, ...
        forceMagnitude, viewStress, viewStrain, defHBool, defScale, defAFace, defFFace, defFType, defFMag)
    
end
%}
