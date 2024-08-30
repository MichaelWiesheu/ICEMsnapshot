
clc
clear all
close all

restoredefaultpath
addpath(genpath("../"))

warning('off','nrbderiv:SecondDerivative')

%% Create geometry
opts.draw_geometry = false;

[srfRotor, patchesRotor, opts] = PMSM_rotor1(opts);
[srfStator, patchesStator, opts] = PMSM_stator1(opts);

%% Initialization
MotorOpt = MotorOptimizationFull();
MotorOpt.Options = opts;
MotorOpt.MachineProperties.Length = 0.1;
MotorOpt.FormFunctionDegree = [2,2];
MotorOpt.DegreeQuadrature = MotorOpt.FormFunctionDegree + 1;
MotorOpt.SubdivisionsPatchesRotor = [3,3];
MotorOpt.SubdivisionsPatchesStator = [3,3];

MotorOpt.importSurface(srfRotor, srfStator);
%% Materials
MatAir = Mat_Air();
MatCopper = Mat_Copper();
MatIron = M_27();
MatMagnet = NdFeB_Br14();

MotorOpt.resetMaterials();
MotorOpt.setMaterial(MatAir, patchesRotor.Air, patchesStator.Air);
MotorOpt.setMaterial(MatCopper, [], patchesStator.Windings);
MotorOpt.setMaterial(MatIron, patchesRotor.Iron, patchesStator.Iron);
MotorOpt.setMaterial(MatMagnet, patchesRotor.Magnets, []);
%%

% %% Boundary conditions
MotorOpt.setDirichletBoundaries([4], [8,12,16,20,24,32])
MotorOpt.setCouplingBoundaries([7,8,9,10,12], [2,3,4,9,10,11,13,14,15,17,18,19,21,22,23,25,26,28]);
MotorOpt.setAntiPeriodicBoundaries([5,2,11], [3,1,6], [27,29,30,31], [1,5,6,7]);
%% Magnet remanence
MotorOpt.resetRemanenceRotor();
MotorOpt.addRemanenceRotor([8], MatMagnet, pi/2);

%%
MotorOpt.XLimits = [-35, 35] * 1e-3;
MotorOpt.YLimits = [10, 72] * 1e-3;
MotorOpt.plotGeometry()
%% Stiffness matrices
MotorOpt.generateStiffnessMatrices();
%% Coupling Spaces and Meshes
MotorOpt.generateCouplingSpacesAndMeshes();
%% Coupling settings
nharm = 15;
mult = 6;
CosValues = 3:mult:mult*nharm/2;
SinValues = 3:mult:mult*nharm/2;

MotorOpt.setCouplingHarmonics(SinValues, CosValues);
MotorOpt.generateCouplingMatrices();

%% Current settings
MotorOpt.resetPhases();
MotorOpt.addPhase("U", [8,10,13,14], 1)
MotorOpt.addPhase("U", [7,9,11,12], 1)
MotorOpt.addPhase("U", [26,28,30,31], 1)
MotorOpt.addPhase("U", [103,105,108,109], -1)

MotorOpt.addPhase("W", [27,29,32,33], -1)
MotorOpt.addPhase("W", [46,48,51,52], -1)
MotorOpt.addPhase("W", [45,47,49,50], -1)
MotorOpt.addPhase("W", [64,66,68,69], -1)

MotorOpt.addPhase("V", [65,67,70,71], 1)
MotorOpt.addPhase("V", [83,85,87,88], 1)
MotorOpt.addPhase("V", [84,86,89,90], 1)
MotorOpt.addPhase("V", [102,104,106,107], 1)

PolePairs  = 3;
NumberWindings = 12;
phase0 = 0;
Iapp = 10;
MotorOpt.setCurrentOptions(Iapp, PolePairs, NumberWindings, phase0);

%%
MotorOpt.setRotorGeometryFunction(@(x) PMSM_rotor1(x));
MotorOpt.setStatorGeometryFunction(@(x) PMSM_stator1(x));
warning('off','all')

MotorOpt.resetOptimizationParameters();

MotorOpt.addOptimizationParameter("ScaleR", 1, 0.5, 2);
MotorOpt.addOptimizationParameter("LENGTH", 0.1, 0.08, 0.12);
MotorOpt.addOptimizationParameter("OPERATING_ANGLE", 0, -20, 20);

MotorOpt.addOptimizationParameter("MAG", 7e-3, 6e-3, 15e-3);
MotorOpt.addOptimizationParameter("MH", 7e-3, 2e-3, 12e-3);
MotorOpt.addOptimizationParameter("MW", 19e-3, 10e-3, 25e-3);

MotorOpt.addOptimizationParameter("SW1", 4e-3, 2e-3, 6e-3);
MotorOpt.addOptimizationParameter("SW2", 2.3e-3, 1e-3, 4e-3);
MotorOpt.addOptimizationParameter("SW3", 1e-3, 0.5e-3, 1.5e-3);
MotorOpt.addOptimizationParameter("SW4", 8.25e-3, 2e-3, 20e-3);
MotorOpt.addOptimizationParameter("SD1", 135e-3, 90e-3, 160e-3);
MotorOpt.addOptimizationParameter("Sr1", 1e-3, 0.5e-3, 2e-3);

MotorOpt.setConstraintFunction(@(x)PMSM_constraints1(x))

MotorOpt.resetOptimizationControlPoints();

MotorOpt.addOptimizationControlPointRotor([43, 83], -0.e-3, -1.5e-3, 0)
MotorOpt.addOptimizationControlPointRotor([44, 84], -0.e-3, -1.5e-3, 0)
MotorOpt.addOptimizationControlPointRotor([40,47,48, 79, 80, 50], -0.e-3, -1e-3, 0) 
MotorOpt.addOptimizationControlPointRotor([59, 66], -0.e-3, -1.5e-3, 0)
MotorOpt.addOptimizationControlPointRotor([60, 65], -0.e-3, -1.5e-3, 0)
MotorOpt.addOptimizationControlPointRotor([61, 64], -0.e-3, -1.5e-3, 0)
MotorOpt.addOptimizationControlPointRotor([62, 63], -0.e-3, -1.5e-3, 0)

MotorOpt.addOptimizationControlPointStator([194, 198, 199, 273,276,277, 348,351,352,423,426,427, 498,501,502,573,576,577], 0, 0, 0.5e-3)
MotorOpt.addOptimizationControlPointStator([196, 202, 274, 280, 349,355,424,430,499,505,574,578], 0, 0, 0.5e-3)
MotorOpt.addOptimizationControlPointStator([195, 201, 279,354,429,504,579], 0, 0, 0.5e-3)



MotorOpt.OptimizationAngles = [0:2:19];
MotorOpt.MinOptiTorque = 1.5;
MotorOpt.OptimizationType = "nonlinear";
MotorOpt.OptimizationMethod = "interior-point";
MotorOpt.MultCost = 0.05;
MotorOpt.MultTsig = 10;
MotorOpt.MultPower = 4;
%%
MotorOpt.setRotationAngle(0);
MotorOpt.setCurrent(Iapp, 0);
MotorOpt.solveMagneticPotentialNewton();
MotorOpt.plotBResulting();
MotorOpt.calcTorqueBrBtRotor()
MotorOpt.calcTorqueBrBtStator()
MotorOpt.calcTorqueCoupling()
%%
xopt = MotorOpt.optimize();
