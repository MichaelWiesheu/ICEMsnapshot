clc
clear all
close all

restoredefaultpath
addpath(genpath("../"))

warning('off','nrbderiv:SecondDerivative')

%%
Date = 'final';
Iter = 120;

OptimizationHistory = load(['plots/' Date '/OptimizationHistory.mat']).OptimizationHistory;
IndIter = find(vertcat(OptimizationHistory.iteration) == Iter);
OptimizationValues = load(['plots/' Date '/OptimizationHistory.mat']).OptimizationHistory(IndIter);

scaleR = OptimizationValues.Parameters(1).OptVal;
length = OptimizationValues.Parameters(2).OptVal;
operatingAngle = OptimizationValues.Parameters(3).OptVal;
currents = OptimizationValues.Currents*scaleR;

%% Create geometry
opts.draw_geometry = false;

[srfRotor, patchesRotor, opts] = PMSM_rotor1(opts);
[srfStator, patchesStator, opts] = PMSM_stator1(opts);

patchesRotor.Iron = [patchesRotor.Iron, 2,5]; % Add iron bridges
patchesRotor.Air = setdiff(patchesRotor.Air, [2,5]); % Remove Air

geoRotor = mp_geo_load(['plots/' Date '/RotorOpt' num2str(Iter) '.txt']);
for iSrf = 1:numel(geoRotor)
    geoRotor(iSrf).nurbs.coefs(1:2,:,:) = geoRotor(iSrf).nurbs.coefs(1:2,:,:)*scaleR;
end
geoStator = mp_geo_load(['plots/' Date '/StatorOpt' num2str(Iter) '.txt']);
for iSrf = 1:numel(geoStator)
    geoStator(iSrf).nurbs.coefs(1:2,:,:) = geoStator(iSrf).nurbs.coefs(1:2,:,:)*scaleR;
end
srfRotor = vertcat(geoRotor.nurbs);
srfStator = vertcat(geoStator.nurbs);


%% Initialization
PMSMmotor = MotorSimulation();
PMSMmotor.MachineProperties.Length = length;
PMSMmotor.FormFunctionDegree = [2,2];
PMSMmotor.DegreeQuadrature = PMSMmotor.FormFunctionDegree + 1;
PMSMmotor.SubdivisionsPatchesRotor = [3,3];
PMSMmotor.SubdivisionsPatchesStator = [3,3];

PMSMmotor.importSurface(srfRotor, srfStator);
%% Materials
MatAir = Mat_Air();
MatCopper = Mat_Copper();
MatIron = M_27();
MatMagnet = NdFeB_Br14();

PMSMmotor.resetMaterials();
PMSMmotor.setMaterial(MatAir, patchesRotor.Air, patchesStator.Air);
PMSMmotor.setMaterial(MatCopper, [], patchesStator.Windings);
PMSMmotor.setMaterial(MatIron, patchesRotor.Iron, patchesStator.Iron);
PMSMmotor.setMaterial(MatMagnet, patchesRotor.Magnets, []);
%%

% %% Boundary conditions
PMSMmotor.setDirichletBoundaries([4], [8,12,16,20,24,32])
PMSMmotor.setCouplingBoundaries([7,8,9,10,12], [2,3,4,9,10,11,13,14,15,17,18,19,21,22,23,25,26,28]);
PMSMmotor.setAntiPeriodicBoundaries([5,2,11], [3,1,6], [27,29,30,31], [1,5,6,7]);
% %% Magnet remanence
PMSMmotor.addRemanenceRotor([8], MatMagnet, pi/2);

%% Stiffness matrices
PMSMmotor.generateStiffnessMatrices();
%% Coupling Spaces and Meshes
PMSMmotor.generateCouplingSpacesAndMeshes();
%% Coupling settings
nharm = 15;
mult = 6;
CosValues = 3:mult:mult*nharm/2;
SinValues = 3:mult:mult*nharm/2;

PMSMmotor.setCouplingHarmonics(SinValues, CosValues);
PMSMmotor.generateCouplingMatrices();

%% Current settings
PMSMmotor.resetPhases();
PMSMmotor.addPhase("U", [8,10,13,14], 1)
PMSMmotor.addPhase("U", [7,9,11,12], 1)
PMSMmotor.addPhase("U", [26,28,30,31], 1)
PMSMmotor.addPhase("U", [103,105,108,109], -1)

PMSMmotor.addPhase("W", [27,29,32,33], -1)
PMSMmotor.addPhase("W", [46,48,51,52], -1)
PMSMmotor.addPhase("W", [45,47,49,50], -1)
PMSMmotor.addPhase("W", [64,66,68,69], -1)

PMSMmotor.addPhase("V", [65,67,70,71], 1)
PMSMmotor.addPhase("V", [83,85,87,88], 1)
PMSMmotor.addPhase("V", [84,86,89,90], 1)
PMSMmotor.addPhase("V", [102,104,106,107], 1)

PolePairs  = 3;
NumberWindings = 12;
phase0 = operatingAngle;
ApplicationCurrent = currents(end);
PMSMmotor.setCurrentOptions(ApplicationCurrent, PolePairs, NumberWindings, phase0);

%%
PMSMmotor.setRotationAngle(0);
PMSMmotor.plotGeometry()
ControlPointsRotor = [43, 83, 44, 84, 40,47,48, 79, 80, 50, 59, 66, 60, 65, 61, 64, 62, 63];
ControlPointsStator = [194, 198, 199, 273,276,277, 348,351,352,423,426,427, 498,501,502,573,576,577,...
    196, 202, 274, 280, 349,355,424,430,499,505,574,578, ...
    195, 201, 279,354,429,504,579];
PMSMmotor.plotControlPoints([1,2], false, [], [], 5, TUDa_getColor("6b"));
PMSMmotor.plotControlPoints([1,2], false, ControlPointsRotor, ControlPointsStator, 25, TUDa_getColor("9b"));

set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
xlim([-1,1]*0.035)
ylim([0.012,0.069])
xticks([-0.03:0.01:0.03])
yticks([0.01:0.01:0.06])
set(gca, "FontSize", 16)
exportgraphics(gca,'GeometryOpt.pdf') 
exportgraphics(gca,'GeometryOpt.png') 

%%
PMSMmotor.setRotationAngle(0);
PMSMmotor.setCurrent(ApplicationCurrent,0);
PMSMmotor.solveMagneticPotentialNewton();
PMSMmotor.plotBResulting()

set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
xlim([-1,1]*0.035)
ylim([0.012,0.069])
xticks([-0.03:0.01:0.03])
yticks([0.01:0.01:0.06])
set(gca, "FontSize", 16)
exportgraphics(gca,'PMSMfieldOpt.pdf','Resolution',1000) 
exportgraphics(gca,'PMSMfieldOpt.png') 

%%

J = currents(end)*PMSMmotor.NumWindings/PMSMmotor.Phases(1).Area

SymMult = 1;
% angles = 0:2:19;
angles = 0:1:19;
Torques = zeros(numel(currents), numel(angles));

for iCurrents = 1:numel(currents)
    for iAngle = 1:numel(angles)
        disp(num2str(angles(iAngle)));
        PMSMmotor.setRotationAngle(angles(iAngle));
        PMSMmotor.setCurrent(currents(iCurrents), angles(iAngle));
        PMSMmotor.solveMagneticPotentialNewton();
        T = SymMult * PMSMmotor.calcTorqueCoupling();
        Torques(iCurrents, iAngle) = T;
    end
end

save("OptTorques.mat", "angles", "Torques", "currents");

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
figure(32)
% clf
hold on
plot(angles, Torques(1,:), "Color", TUDa_getColor("8d"), LineWidth=2);
plot(angles, Torques(2,:), "Color", TUDa_getColor("7a"));
plot(angles, Torques(3,:), "Color", TUDa_getColor("9a"));
% ylim([0,1.9])
% legend("Coupling Rt", "BrBt Rotor", "BrBt Stator")
grid on
xlabel("Rotation angle")
ylabel("Torque (Nm/m)")
set(gca, "FontSize",12)

%%

for iMat = 1:numel(PMSMmotor.Materials)
    patchesRt = PMSMmotor.Materials(iMat).PatchesRotor;
    patchesSt = PMSMmotor.Materials(iMat).PatchesStator;
    rho = PMSMmotor.Materials(iMat).Material.getRho();
    cost = PMSMmotor.Materials(iMat).Material.getCost();
    Costi(iMat) = 0;
    Costi(iMat) = Costi(iMat) +  rho * cost * length * sum(PMSMmotor.PatchAreaRotor(patchesRt));
    Costi(iMat) = Costi(iMat) +  rho * cost * length * sum(PMSMmotor.PatchAreaStator(patchesSt));
end
Cost = sum(Costi);


Tmean = mean(Torques, 2)
Tstd = std(Torques, 1, 2)

sigma = 59600000;
PowerLoss = J^2*PMSMmotor.Phases(1).Area/PMSMmotor.NumWindings/sigma*length;

m1 = OptimizationValues.Mult.MultCost;
m2 = OptimizationValues.Mult.MultTsig;
m3 = OptimizationValues.Mult.MultPower;

fopt = m1*Cost + m2*sum(Tstd) + m3*PowerLoss;


%%
currentsMapOpt = linspace(0, max(currents), 21);
PhaseAnglesOpt = -180:5:180;%linspace(0,90,20);
TorquesMapOpt = zeros(numel(currents), numel(PhaseAnglesOpt), numel(angles));

for iCurrents = 1:numel(currentsMapOpt)
    for iPhase = 1:numel(PhaseAnglesOpt)
        PMSMmotor.setCurrentOptions(ApplicationCurrent, PolePairs, NumberWindings, PhaseAnglesOpt(iPhase));
        for iAngle = 1:numel(angles)
%             disp(num2str(angles(iAngle)));
            PMSMmotor.setRotationAngle(angles(iAngle));
            PMSMmotor.setCurrent(currentsMapOpt(iCurrents), angles(iAngle));
            PMSMmotor.solveMagneticPotentialNewton();
            T = SymMult * PMSMmotor.calcTorqueCoupling();
            TorquesMapOpt(iCurrents, iPhase, iAngle) = T;
        end
    end
    iCurrents
end

save("OptValues.mat", "TorquesMapOpt", "currentsMapOpt", "PhaseAnglesOpt", "angles");






