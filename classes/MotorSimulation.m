classdef MotorSimulation < handle

    properties
        % Quality of computation
        FormFunctionDegree = [2, 2];                                    % Degree of form functions 1-linear, 2-quadratic
        DegreeQuadrature = [3, 3];                                      % Degree of quadrature, usually use FormFunctionDegree
        SubdivisionsPatchesRotor = [10, 10];                            % Subdivision for each patch, the more the accurate
        SubdivisionsPatchesStator = [10, 10];  

        % Domain Definitions for Rotor and Stator
        SurfaceRotor, SurfaceStator;                                    % NURBS representation of single patches
        GeometryRotor, GeometryStator;                                  % Parametrizatoin of the mapping to the physical domain and derivatives
        BoundariesRotor, BoundariesStator                               % Information about the patches with an outer boundary
        InterfacesRotor, InterfacesStator                               % Information about the connections between patches
        BoundaryInterfacesRotor, BoundaryInterfacesStator;              % Boundary Information
        MultipatchSpaceRotor, MultipatchSpaceStator;                    % Space stores the shape functions of discrete spaces 
        MultipatchSpaceRotorGeo, MultipatchSpaceStatorGeo;   
        MultipatchMeshRotor, MultipatchMeshStator                       % Mesh stores partition of domain and quadrature rule
        MpSpaceEvalRt, MpSpaceEvalSt;                                   % Evaluated Structure
        MpSpaceGeoEvalRt, MpSpaceGeoEvalSt;                             %
        MpMeshEvalRt, MpMeshEvalSt;                                     %

        % Coupling definitions, analogous to rotor and stator for coupling domain
        SurfaceRotorCoupling, SurfaceStatorCoupling;
        GeometryRotorCoupling, GeometryStatorCoupling;
        BoundariesCouplingRotor, BoundariesCouplingStator;              % Coupling Boundaries
        BoundariesDirichletRotor, BoundariesDirichletStator;            % Dirichlet Boundaries   
        BoundariesAntiperiodicRotor, BoundariesAntiperiodicStator;      % Antiperiodic Boundaries 
        InterfacesRotorCoupling, InterfacesStatorCoupling;
        BoundaryInterfacesRotorCoupling, BoundaryInterfacesStatorCoupling;
        MultipatchSpaceRotorCoupling, MultipatchSpaceStatorCoupling;    
        MultipatchMeshRotorCoupling, MultipatchMeshStatorCoupling;
        BoundariesRotorCoupling, BoundariesStatorCoupling;
        
        % DOFS;
        NumberPatchesStator, NumberPatchesRotor;                        % Patch counter
        NumberDOFRotor, NumberDOFStator;                                % Number of degrees of freedom
        NumberDOFRotorCoupling, NumberDOFStatorCoupling;                % Number DOF for coupling
        DirichletDOFSRotor, DirichletDOFSStator;                        % Dirichlet degrees of freedom
        IndependentDOFsRotor, IndependentDOFsStator;                    % Independent Dofs (not dirichlet)
        APdofsRotorLeft, APdofsRotorRight;                              % Antiperiodic DOFS Rotor
        APdofsStatorLeft, APdofsStatorRight;                            % Antiperiodic DOFS Stator
        CouplingDOFSRotor, CouplingDOFSStator;                          % Coupling Dofs
        MagneticPotentialRotor, MagneticPotentialStator;                % Multiplier for shape function defining the resulting magnetic potential A from which real B-values are derived

        % Stiffnes Matrices
        StiffMatsRotor, StiffMatsStator;                                % Cell array containting stiffness matrices for different materials
        StiffMatRotor, StiffMatStator;                                  % Assembled Stiffness Matrices for rotor and stator
        MassMatsRotor, MassMatsStator;
        MassMatRotor, MassMatStator;
        RHSRotor, RHSStator;                                            % The right hand side of equations, containing the permanent magnet
        dRHSRotor, dRHSStator;                                          % Derivative of rhs after Iapp
        dRHSStatorangle, dRHSRotorangle                                 % derivative after Angle
        ddRHSStatorangle, ddRHSRotorrangle                               %secobd derivative after angle
        % Coupling Matrices
        NumberDOFCoupling;
        HarmonicsCos, HarmonicsSin, HarmonicsAll, CouplingIndices;
        RotationAngle = 0;
        CouplingMatrixRotor, CouplingMatrixStator;
        CouplingMatrixRotorInit, CouplingMatrixStatorInit;
        CouplingMatrixRotorFull, CouplingMatrixStatorFull;
%         CouplingMatrixStatorFullSin, CouplingMatrixStatorFullCos;
        RotationMatrix, RotationMatrixDer, RotationMatrixDer2, RotationMatrixDerFull;                              % Rotation matrix and derivative wrt rotation angle
        CouplingValues;

        CouplingMeshesRotor, CouplingMeshesStator;                      % Only for
        CouplingSpacesRotor, CouplingSpacesStator;                      % Torque
        CouplingSpacesRotorEval, CouplingSpacesStatorEval;              % Evaluation so far!

        % Electric currents
        ApplicationCurrent, PolePairs, NumWindings, Phase0;
        Phases;

        % Material properties
        Materials;
        Remanence;
        MuVac = 4*pi*1e-7;
        MachineProperties;                                              % Struct with length and material weights
        PatchAreaRotor, PatchAreaStator;                                % Area of each patch
        PatchInertiaRotor, PatchInertiaStator;
        NewtonResidual = 1e-6;

        % Forces
        ForceDensityRotor, ForceDensityStator;
        
        % Plot definitions
        ParamPlotPointsQuiverRotor, ParamPlotPointsQuiverStator;        % Parametric points for plotting B-Field arrows
        GeoPlotPointsQuiverRotor, GeoPlotPointsQuiverStator;            % Physical points for plotting B-Field arrows
        RotorShapeValuesQuiver, StatorShapeValuesQuiver;                % Shape values and gradients for arrow plot points
        ParamPlotPointsRotor, ParamPlotPointsStator;                    % Parametric points for surface evaluation
        GeoPlotPointsRotor, GeoPlotPointsStator;                        % Physical points for surface evaluation
        RotorShapeValues, StatorShapeValues;                            % Shape values and gradients for surface plot points
        NumberQuiversRotor = 300, NumberQuiversStator = 800;            % Number of Arrows for drawing the magnet field
        PlotRemanence = {};                                             % Struct for remanent patches
        firstBplot=true;                                                % indicate if plot points must be calculated

        % Values for Evaluation on single points
        PlotResolutionX = 4;
        PlotResolutionY = 4;
        EvalPointsRotor, EvalPointsStator;                        % Parametric points at which the Air gap values are evaluated
        EvalPatchesRotor, EvalPatchesStator;                      % Corresponding patch for each point in GapEvalPoints
        EvalShapeRotor, EvalShapeStator;                            % Shape values and gradients at evaluation points

        % Handles for Figure Geometry
        FigureGeometry;                                                     % Figure handle
        FigureGeometryRotor, FigureGeometryStator;          
        FigureBoundariesRotor, FigureBoundariesStator;
        FigureInterfacesRotor, FigureInterfacesStator;
        % Handles for Figure Potential Lines
        FigurePotLines;                                                     % Figure handle
        FigurePotLinesGeometryRotor, FigurePotLinesGeometryStator;          % Surf handle
        FigurePotLinesBoundariesRotor, FigurePotLinesBoundariesStator;
        FigurePotLinesInterfacesRotor, FigurePotLinesInterfacesStator;
        FigurePotLinesPotentialLinesRotor, FigurePotLinesPotentialLinesStator;
        FigurePotLinesMin, FigurePotLinesMax;
        % Handles for Figure of B Resulting
        FigureBRes;
        FigureBResGeometryRotor, FigureBResGeometryStator;
        FigureBResBoundariesRotor, FigureBResBoundariesStator;
        FigureBResInterfacesRotor, FigureBResInterfacesStator;
        FigureBResBMagRotor, FigureBResBMagStator;
        FigureBResQuiver;
        %
        QuiverPlot;                                     % temp value for storing the arrows of remanence
        XLimits = [-Inf, Inf];                      % limits for Plots
        YLimits = [-Inf, Inf];                      % limits for plots
        PlotColor;                                      % Struct with colors for Air, Iron, Magnet
        PlotFaceAlpha = 0.8;
        NumPotLines = 40;                               % Number of potential lines for potential line plot

        ControlPointsRotor, ControlPointsStator;
        ControlPointsRotorOffset;
    end

    methods
        %% constructor, set material properties
        function obj = MotorSimulation()

            obj.MachineProperties.Length = 1;

            obj.resetMaterials();
            obj.resetPhases();
            
            obj.PlotColor.Lines = [0.5, 0.5, 0.5];
            obj.PlotColor.Contour = [0.2, 0.2, 0.2];

        end

        %% Import the Surface, create Geometry, Meshes and Spaces
        function importSurface(obj, surfaceRotor, surfaceStator)  
            warning('off', 'nrbderiv:SecondDerivative')        
            % Rotor
            disp("Creating Rotor spaces and meshes...")
            obj.SurfaceRotor = surfaceRotor;            
            [obj.GeometryRotor, obj.BoundariesRotor, obj.InterfacesRotor, ~, obj.BoundaryInterfacesRotor] = mp_geo_load(obj.SurfaceRotor, 1e-8);
            obj.NumberPatchesRotor = numel(obj.GeometryRotor);
            meshes = cell (1, obj.NumberPatchesRotor);
            spaces  = cell (1, obj.NumberPatchesRotor);
            spaces_geo = cell (1, obj.NumberPatchesRotor);
            % create space partition and form functions
            for iptc = 1:obj.NumberPatchesRotor
              [knots, zeta] = kntrefine (obj.GeometryRotor(iptc).nurbs.knots, obj.SubdivisionsPatchesRotor-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
              rule      = msh_gauss_nodes (obj.DegreeQuadrature);
              [qn, qw]  = msh_set_quad_nodes (zeta, rule);
              meshes{iptc} = msh_cartesian (zeta, qn, qw, obj.GeometryRotor(iptc));
              spaces{iptc}  = sp_bspline (knots, obj.FormFunctionDegree, meshes{iptc});
              spaces_geo{iptc} = sp_nurbs(obj.GeometryRotor(iptc).nurbs, meshes{iptc});
            end
            obj.MultipatchMeshRotor = msh_multipatch(meshes, obj.BoundariesRotor);
            obj.MultipatchSpaceRotor  = sp_multipatch(spaces, obj.MultipatchMeshRotor, obj.InterfacesRotor, obj.BoundaryInterfacesRotor);
            obj.MultipatchSpaceRotorGeo = sp_multipatch(spaces_geo, obj.MultipatchMeshRotor, obj.InterfacesRotor, obj.BoundaryInterfacesRotor);
            
            % Pre-Evaluate Spaces and Meshes
            obj.MpMeshEvalRt = cell (1, obj.NumberPatchesRotor);
            obj.MpSpaceEvalRt  = cell (1, obj.NumberPatchesRotor);
            obj.MpSpaceGeoEvalRt = cell (1, obj.NumberPatchesRotor);
            for iPatch = 1:obj.NumberPatchesRotor
                obj.MpMeshEvalRt{iPatch} = msh_precompute(obj.MultipatchMeshRotor.msh_patch{iPatch});
                obj.MpSpaceEvalRt{iPatch}  = sp_precompute(obj.MultipatchSpaceRotor.sp_patch{iPatch}, obj.MpMeshEvalRt{iPatch}, 'value', true, 'gradient', true);
                obj.MpSpaceGeoEvalRt{iPatch}  = sp_precompute_param(obj.MultipatchSpaceRotorGeo.sp_patch{iPatch}, obj.MpMeshEvalRt{iPatch}, 'value', true, 'gradient', true);
            end
            % init values
            obj.NumberDOFRotor = obj.MultipatchSpaceRotor.ndof;
            obj.RHSRotor = zeros(obj.NumberDOFRotor, 1);
            obj.MagneticPotentialRotor = zeros(obj.NumberDOFRotor, 1);
            % Rotor Control Points
            obj.ControlPointsRotor = zeros(obj.MultipatchSpaceRotorGeo.ndof, 4);
            for iPatch  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{iPatch};
                obj.ControlPointsRotor(ind_loc, 1) = reshape(obj.GeometryRotor(iPatch).nurbs.coefs(1, :, :, :)./obj.GeometryRotor(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPointsRotor(ind_loc, 2) = reshape(obj.GeometryRotor(iPatch).nurbs.coefs(2, :, :, :)./obj.GeometryRotor(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPointsRotor(ind_loc, 3) = reshape(obj.GeometryRotor(iPatch).nurbs.coefs(3, :, :, :)./obj.GeometryRotor(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPointsRotor(ind_loc, 4) = reshape(obj.GeometryRotor(iPatch).nurbs.coefs(4, :, :, :), [], 1);
            end
            obj.ControlPointsRotorOffset = zeros(size(obj.ControlPointsRotor));

            % Stator
            disp("Creating Stator spaces and meshes...")
            obj.SurfaceStator = surfaceStator;
            [obj.GeometryStator, obj.BoundariesStator, obj.InterfacesStator, ~, obj.BoundaryInterfacesStator] = mp_geo_load(obj.SurfaceStator);
            obj.NumberPatchesStator = numel(obj.GeometryStator);
            meshes = cell (1, obj.NumberPatchesStator);
            spaces  = cell (1, obj.NumberPatchesStator);
            spaces_geo = cell (1, obj.NumberPatchesStator);
            % create space partition and form functions
            for iptc = 1:obj.NumberPatchesStator
              [knots, zeta] = kntrefine (obj.GeometryStator(iptc).nurbs.knots, obj.SubdivisionsPatchesStator-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
              rule      = msh_gauss_nodes (obj.DegreeQuadrature);
              [qn, qw]  = msh_set_quad_nodes (zeta, rule);
              meshes{iptc} = msh_cartesian (zeta, qn, qw, obj.GeometryStator(iptc));
              spaces{iptc}  = sp_bspline (knots, obj.FormFunctionDegree, meshes{iptc});
              spaces_geo{iptc} = sp_nurbs(obj.GeometryStator(iptc).nurbs, meshes{iptc});
            end
            obj.MultipatchMeshStator = msh_multipatch(meshes, obj.BoundariesStator);
            obj.MultipatchSpaceStator  = sp_multipatch(spaces, obj.MultipatchMeshStator, obj.InterfacesStator, obj.BoundaryInterfacesStator);
            obj.MultipatchSpaceStatorGeo = sp_multipatch(spaces_geo, obj.MultipatchMeshStator, obj.InterfacesStator, obj.BoundaryInterfacesStator);
            % Pre-Evaluate Spaces and Meshes
            obj.MpMeshEvalSt = cell (1, obj.NumberPatchesStator);
            obj.MpSpaceEvalSt  = cell (1, obj.NumberPatchesStator);
            obj.MpSpaceGeoEvalSt = cell (1, obj.NumberPatchesStator);
            for iPatch = 1:obj.NumberPatchesStator
                obj.MpMeshEvalSt{iPatch} = msh_precompute(obj.MultipatchMeshStator.msh_patch{iPatch});
                obj.MpSpaceEvalSt{iPatch}  = sp_precompute(obj.MultipatchSpaceStator.sp_patch{iPatch}, obj.MpMeshEvalSt{iPatch}, 'value', true, 'gradient', true);
                obj.MpSpaceGeoEvalSt{iPatch}  = sp_precompute_param(obj.MultipatchSpaceStatorGeo.sp_patch{iPatch}, obj.MpMeshEvalSt{iPatch}, 'value', true, 'gradient', true);
            end
            % init values
            obj.NumberDOFStator = obj.MultipatchSpaceStator.ndof;
            obj.RHSStator = zeros(obj.NumberDOFStator, 1);
            obj.MagneticPotentialStator = zeros(obj.NumberDOFStator, 1);
            % Stator Control Points
            obj.ControlPointsStator = zeros(obj.MultipatchSpaceStatorGeo.ndof, 4);
            for iPatch  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{iPatch};
                obj.ControlPointsStator(ind_loc, 1) = reshape(obj.GeometryStator(iPatch).nurbs.coefs(1, :, :, :)./obj.GeometryStator(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPointsStator(ind_loc, 2) = reshape(obj.GeometryStator(iPatch).nurbs.coefs(2, :, :, :)./obj.GeometryStator(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPointsStator(ind_loc, 3) = reshape(obj.GeometryStator(iPatch).nurbs.coefs(3, :, :, :)./obj.GeometryStator(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPointsStator(ind_loc, 4) = reshape(obj.GeometryStator(iPatch).nurbs.coefs(4, :, :, :), [], 1);
            end

            % Calc properties and plot points
            obj.calculateParamPlotPoints();
            obj.calculateAreaAndInertia();
            obj.calculateQuiverPlotPoints();
            obj.calculateGeoPlotPoints();
        end

        % Material definitions
        function setMaterial(obj, material, patchesRotor, patchesStator)
            pos = numel(obj.Materials)+1;
            obj.Materials(pos).Material = material;
            obj.Materials(pos).PatchesRotor = patchesRotor;
            obj.Materials(pos).PatchesStator = patchesStator;

            % Calculating Inertia and Mass
            obj.Materials(pos).MassRotor = obj.MachineProperties.Length * sum(obj.PatchAreaRotor(obj.Materials(pos).PatchesRotor)) * obj.Materials(pos).Material.getRho();
            obj.Materials(pos).MassStator = obj.MachineProperties.Length * sum(obj.PatchAreaStator(obj.Materials(pos).PatchesStator)) * obj.Materials(pos).Material.getRho();
            obj.Materials(pos).InertiaRotor = obj.MachineProperties.Length * sum(obj.PatchInertiaRotor(obj.Materials(pos).PatchesRotor)) * obj.Materials(pos).Material.getRho();
            obj.Materials(pos).InertiaStator = obj.MachineProperties.Length * sum(obj.PatchInertiaStator(obj.Materials(pos).PatchesStator)) * obj.Materials(pos).Material.getRho();
        end

        function resetMaterials(obj)
            obj.Materials = struct('Material', {}, 'PatchesRotor', {}, 'PatchesStator', {}, ...
                'MassRotor', {}, 'MassStator', {}, 'InertiaRotor', {}, 'InertiaStator', {});
        end

        %% Set the evaluation points in the air gap by specifying parametric coortinate points and their respective patch
        function setEvaluationPointsRotor(obj, evalPoints, evalPatches)
            obj.EvalPointsRotor = evalPoints;
            obj.EvalPatchesRotor = evalPatches;    
            obj.EvalShapeRotor = shp_eval_mp(obj.MultipatchSpaceRotor, obj.GeometryRotor, evalPoints, evalPatches, {'value', 'gradient'});
        end

        function setEvaluationPointsStator(obj, evalPoints, evalPatches)
            obj.EvalPointsStator = evalPoints;
            obj.EvalPatchesStator = evalPatches;    
            obj.EvalShapeStator = shp_eval_mp(obj.MultipatchSpaceStator, obj.GeometryStator, evalPoints, evalPatches, {'value', 'gradient'});
        end

        function [brt, bst] = evalEvalPoints(obj)
            brt = []; bst = [];
            if ~isempty(obj.EvalShapeRotor)
                gradNi_ui = double(ttv(obj.EvalShapeRotor{2}, obj.MagneticPotentialRotor, 3));
                brt = sqrt(sum(gradNi_ui.^2, 2));
            end
            if ~isempty(obj.EvalShapeStator)
                gradNi_ui = double(ttv(obj.EvalShapeStator{2}, obj.MagneticPotentialStator, 3));
                bst = sqrt(sum(gradNi_ui.^2, 2));
            end
        end

        %% Dirichlet Sides, will be set to homogenious Dirichlet conditions
        function setDirichletBoundaries(obj, dirichletSidesRotor, dirichletSidesStator)
            obj.BoundariesDirichletRotor = dirichletSidesRotor;
            obj.BoundariesDirichletStator = dirichletSidesStator;
            obj.DirichletDOFSRotor = [];
            obj.DirichletDOFSStator = [];
            for side = dirichletSidesRotor
                LocalDOFS = obj.MultipatchSpaceRotor.boundary.gnum{side};
                GlobalDOFS = obj.MultipatchSpaceRotor.boundary.dofs(LocalDOFS);
                obj.DirichletDOFSRotor = union(obj.DirichletDOFSRotor, GlobalDOFS);
            end    
            for side = dirichletSidesStator
                LocalDOFS = obj.MultipatchSpaceStator.boundary.gnum{side};
                GlobalDOFS = obj.MultipatchSpaceStator.boundary.dofs(LocalDOFS);
                obj.DirichletDOFSStator = union(obj.DirichletDOFSStator, GlobalDOFS);
            end  
            obj.IndependentDOFsRotor = setdiff([1:obj.NumberDOFRotor]', obj.DirichletDOFSRotor);
            obj.IndependentDOFsStator = setdiff([1:obj.NumberDOFStator]', obj.DirichletDOFSStator);
        end

        %% Coupling Sides, will be coupled for rotation
        function setCouplingBoundaries(obj, couplingSidesRotor, couplingSidesStator)
            obj.BoundariesCouplingRotor = couplingSidesRotor;
            obj.BoundariesCouplingStator = couplingSidesStator;
            obj.CouplingDOFSRotor = [];
            obj.CouplingDOFSStator = [];
            for sides = couplingSidesRotor
                LocalDOFS = obj.MultipatchSpaceRotor.boundary.gnum{sides};
                GlobalDOFS = obj.MultipatchSpaceRotor.boundary.dofs(LocalDOFS);
                obj.CouplingDOFSRotor = union(obj.CouplingDOFSRotor, GlobalDOFS);
            end    
            for sides = couplingSidesStator
                LocalDOFS = obj.MultipatchSpaceStator.boundary.gnum{sides};
                GlobalDOFS = obj.MultipatchSpaceStator.boundary.dofs(LocalDOFS);
                obj.CouplingDOFSStator = union(obj.CouplingDOFSStator, GlobalDOFS);
            end  
        end
        %% Antiperiodic Sides, for rotor and stator
        function setAntiPeriodicBoundaries(obj, APSideRotorLeft, APSideRotorRight, APSideStatorLeft, APSideStatorRight)
            obj.BoundariesAntiperiodicRotor = union(APSideRotorLeft, APSideRotorRight);
            obj.BoundariesAntiperiodicStator = union(APSideStatorLeft, APSideStatorRight);
            % ROTOR
            obj.APdofsRotorLeft = [];
            obj.APdofsRotorRight = [];
            for side = APSideRotorLeft
                LocalDOFS = obj.MultipatchSpaceRotor.boundary.gnum{side};
                GlobalDOFS = obj.MultipatchSpaceRotor.boundary.dofs(LocalDOFS);
                obj.APdofsRotorLeft = union(obj.APdofsRotorLeft, GlobalDOFS);
            end
            obj.APdofsRotorLeft = setdiff(obj.APdofsRotorLeft, obj.CouplingDOFSRotor);
            for side = APSideRotorRight
                LocalDOFS = obj.MultipatchSpaceRotor.boundary.gnum{side};
                GlobalDOFS = obj.MultipatchSpaceRotor.boundary.dofs(LocalDOFS);
                obj.APdofsRotorRight = union(obj.APdofsRotorRight, GlobalDOFS);
            end
            obj.APdofsRotorRight = setdiff(obj.APdofsRotorRight, obj.CouplingDOFSRotor);

            % STATOR
            obj.APdofsStatorLeft = [];
            obj.APdofsStatorRight = [];
            for side = APSideStatorLeft
                LocalDOFS = obj.MultipatchSpaceStator.boundary.gnum{side};
                GlobalDOFS = obj.MultipatchSpaceStator.boundary.dofs(LocalDOFS);
                obj.APdofsStatorLeft = union(obj.APdofsStatorLeft, GlobalDOFS);
            end
            obj.APdofsStatorLeft = setdiff(obj.APdofsStatorLeft, obj.CouplingDOFSStator);
            for side = APSideStatorRight
                LocalDOFS = obj.MultipatchSpaceStator.boundary.gnum{side};
                GlobalDOFS = obj.MultipatchSpaceStator.boundary.dofs(LocalDOFS);
                obj.APdofsStatorRight = union(obj.APdofsStatorRight, GlobalDOFS);
            end
            obj.APdofsStatorRight = setdiff(obj.APdofsStatorRight, obj.CouplingDOFSStator);

        end

        function generateCouplingSpacesAndMeshes(obj)
            % Rotor
            for bnd = 1:numel(obj.BoundariesCouplingRotor)
                patch = obj.BoundariesRotor(obj.BoundariesCouplingRotor(bnd)).patches;
                side = obj.BoundariesRotor(obj.BoundariesCouplingRotor(bnd)).faces;
                obj.CouplingMeshesRotor{bnd} = msh_eval_boundary_side(obj.MultipatchMeshRotor.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(obj.MultipatchMeshRotor.msh_patch{patch}, side);
                obj.CouplingSpacesRotor{bnd} = obj.MultipatchSpaceRotor.sp_patch{patch}.constructor(msh_side_int);
                obj.CouplingSpacesRotorEval{bnd} = sp_precompute(obj.CouplingSpacesRotor{bnd}, msh_side_int, 'value', true, 'gradient', true);
            end
            % Stator
            for bnd = 1:numel(obj.BoundariesCouplingStator)
                patch = obj.BoundariesStator(obj.BoundariesCouplingStator(bnd)).patches;
                side = obj.BoundariesStator(obj.BoundariesCouplingStator(bnd)).faces;
                obj.CouplingMeshesStator{bnd} = msh_eval_boundary_side(obj.MultipatchMeshStator.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(obj.MultipatchMeshStator.msh_patch{patch}, side);
                obj.CouplingSpacesStator{bnd} = obj.MultipatchSpaceStator.sp_patch{patch}.constructor(msh_side_int);
                obj.CouplingSpacesStatorEval{bnd} = sp_precompute(obj.CouplingSpacesStator{bnd}, msh_side_int, 'value', true, 'gradient', true);
            end
        end

        function setCouplingHarmonics(obj, sinValues, cosValues)
            obj.HarmonicsSin = reshape(unique(sinValues), 1, []);
            obj.HarmonicsCos = reshape(unique(cosValues), 1, []);
            obj.HarmonicsAll = 0:max([obj.HarmonicsSin, obj.HarmonicsCos]);
            obj.NumberDOFCoupling = numel(obj.HarmonicsSin) + numel(obj.HarmonicsCos);
            obj.CouplingIndices = union(obj.HarmonicsSin*2 + 1, obj.HarmonicsCos*2 + 2);
        end

        function generateCouplingMatrices(obj)
            disp("Generating Coupling Matrices...")
            % Function definitions values
            SinPhi =  @(x, y) sin(reshape(obj.HarmonicsAll, [], 1).* atan2 (reshape(y, [1, size(y)]), reshape(x, [1, size(x)])));
            CosPhi =  @(x, y) cos(reshape(obj.HarmonicsAll, [], 1).* atan2 (reshape(y, [1, size(y)]), reshape(x, [1, size(x)])));
            
            % Rotor Coupling Values
            B_rt = sparse (obj.NumberDOFRotor, 2*numel(obj.HarmonicsAll));
            for sides = 1:numel(obj.BoundariesCouplingRotor)
                dofsRotor = obj.MultipatchSpaceRotor.gnum{obj.BoundariesRotor(obj.BoundariesCouplingRotor(sides)).patches};
                B_rt(dofsRotor, 1:2:end) = B_rt(dofsRotor, 1:2:end) + op_fs_v (obj.CouplingSpacesRotorEval{sides}, obj.CouplingMeshesRotor{sides}, SinPhi);
                B_rt(dofsRotor, 2:2:end) = B_rt(dofsRotor, 2:2:end) + op_fs_v (obj.CouplingSpacesRotorEval{sides}, obj.CouplingMeshesRotor{sides}, CosPhi);
            end
            
            % Rotor Coupling Values
            B_st = sparse (obj.NumberDOFStator, 2*numel(obj.HarmonicsAll));
            for sides = 1:numel(obj.BoundariesCouplingStator)
                dofsStator = obj.MultipatchSpaceStator.gnum{obj.BoundariesStator(obj.BoundariesCouplingStator(sides)).patches};
                B_st(dofsStator, 1:2:end) = B_st(dofsStator, 1:2:end) + op_fs_v (obj.CouplingSpacesStatorEval{sides}, obj.CouplingMeshesStator{sides}, SinPhi);
                B_st(dofsStator, 2:2:end) = B_st(dofsStator, 2:2:end) + op_fs_v (obj.CouplingSpacesStatorEval{sides}, obj.CouplingMeshesStator{sides}, CosPhi);
            end
            
            B_rt = B_rt / obj.MuVac;
            B_st = B_st / obj.MuVac;

            obj.CouplingMatrixRotorFull = B_rt;
            obj.CouplingMatrixStatorFull = B_st;
            obj.CouplingValues = zeros(obj.NumberDOFCoupling, 1);

            obj.CouplingMatrixRotorInit = obj.CouplingMatrixRotorFull(:, obj.CouplingIndices);
            obj.CouplingMatrixStatorInit = obj.CouplingMatrixStatorFull(:, obj.CouplingIndices);

            obj.CouplingMatrixRotor = obj.CouplingMatrixRotorInit;
            obj.CouplingMatrixStator = obj.CouplingMatrixStatorInit;

            obj.setRotationAngle(obj.RotationAngle);
        end

        function setRotationAngle(obj, angle)
            obj.RotationAngle = deg2rad(angle);
            R = sparse(2*numel(obj.HarmonicsAll), 2*numel(obj.HarmonicsAll));
            Rder = sparse(2*numel(obj.HarmonicsAll), 2*numel(obj.HarmonicsAll));
            Rder2 = sparse(2*numel(obj.HarmonicsAll), 2*numel(obj.HarmonicsAll));
            for imode = obj.HarmonicsAll
                Harmonic = obj.HarmonicsAll(imode+1);
                R(2*imode+1, 2*imode+1) =  cos(Harmonic* obj.RotationAngle);
                R(2*imode+2, 2*imode+1) = -sin(Harmonic* obj.RotationAngle);
                R(2*imode+1, 2*imode+2) =  sin(Harmonic* obj.RotationAngle);
                R(2*imode+2, 2*imode+2) =  cos(Harmonic* obj.RotationAngle);

                Rder(2*imode+1, 2*imode+1) = -Harmonic*sin(Harmonic* obj.RotationAngle);
                Rder(2*imode+2, 2*imode+1) = -Harmonic*cos(Harmonic* obj.RotationAngle);
                Rder(2*imode+1, 2*imode+2) =  Harmonic*cos(Harmonic* obj.RotationAngle);
                Rder(2*imode+2, 2*imode+2) = -Harmonic*sin(Harmonic* obj.RotationAngle);

                Rder2(2*imode+1, 2*imode+1) = -Harmonic^2*cos(Harmonic* obj.RotationAngle);
                Rder2(2*imode+2, 2*imode+1) =  Harmonic^2*sin(Harmonic* obj.RotationAngle);
                Rder2(2*imode+1, 2*imode+2) = -Harmonic^2*sin(Harmonic* obj.RotationAngle);
                Rder2(2*imode+2, 2*imode+2) = -Harmonic^2*cos(Harmonic* obj.RotationAngle);
            end
            
            obj.RotationMatrix = R(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDer = Rder(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDer2 = Rder2(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDerFull = Rder;

            Bst = obj.CouplingMatrixStatorFull*R;

            obj.CouplingMatrixStator = Bst(:, obj.CouplingIndices);
        end

        %% Stiffness linear Matrices, with constant mu_r
        function generateStiffnessMatrices(obj)
            disp("Generating Stiffness Matrices...")
            obj.StiffMatRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
            obj.StiffMatStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
            for iMat = 1:numel(obj.Materials)
                % here always linear calculation
                nu = obj.Materials(iMat).Material.getNuLinear();
                obj.StiffMatsRotor{iMat} = op_gradu_gradv_mp_eval(obj.MultipatchSpaceRotor, obj.MpSpaceEvalRt, obj.MultipatchSpaceRotor, obj.MpSpaceEvalRt, obj.MultipatchMeshRotor, obj.MpMeshEvalRt, nu, obj.Materials(iMat).PatchesRotor);
                obj.StiffMatsStator{iMat} = op_gradu_gradv_mp_eval(obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MultipatchMeshStator, obj.MpMeshEvalSt, nu, obj.Materials(iMat).PatchesStator);
                obj.StiffMatRotor = obj.StiffMatRotor + obj.StiffMatsRotor{iMat};
                obj.StiffMatStator = obj.StiffMatStator + obj.StiffMatsStator{iMat};
            end
        end

        function generateMassMatrices(obj)
            disp("Generating Mass Matrices...")
            obj.MassMatRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
            obj.MassMatStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
            for iMat = 1:numel(obj.Materials)
                % here always linear calculation
                sigma = obj.Materials(iMat).Material.getSigma();
                obj.MassMatsRotor{iMat} = op_u_v_mp_eval(obj.MultipatchSpaceRotor, obj.MpSpaceEvalRt, obj.MultipatchSpaceRotor, obj.MpSpaceEvalRt, obj.MultipatchMeshRotor, obj.MpMeshEvalRt, sigma, obj.Materials(iMat).PatchesRotor);
                obj.MassMatsStator{iMat} = op_u_v_mp_eval(obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MultipatchMeshStator, obj.MpMeshEvalSt, sigma, obj.Materials(iMat).PatchesStator);
                obj.MassMatRotor = obj.MassMatRotor + obj.MassMatsRotor{iMat};
                obj.MassMatStator = obj.MassMatStator + obj.MassMatsStator{iMat};
            end
        end

        function resetRemanenceRotor(obj)
            obj.RHSRotor = zeros(size(obj.RHSRotor));
            obj.PlotRemanence = {};
            obj.Remanence = struct('Patch', {}, 'Material', {}, 'InitAngle', {}, 'Relation', {}, 'Description', {});
        end

        %% add permanent magnet properties inside domain
        function addRemanenceRotor(obj, patch, material, angle, relation, description)
            pos = numel(obj.Remanence) + 1;
            obj.Remanence(pos).Patch = patch;
            obj.Remanence(pos).Material = material;
            obj.Remanence(pos).InitAngle = angle;
            if nargin < 5
                relation = angle; description = [];
            end
            obj.Remanence(pos).Relation = relation;
            obj.Remanence(pos).Description = description;

            Hc = material.Br / material.getMuLinear();

            obj.RHSRotor = obj.RHSRotor + Hc*op_gradv_n_bot_mp(obj.MultipatchSpaceRotor, obj.MultipatchMeshRotor, angle, patch);
            % for geomerty plot
            obj.PlotRemanence{patch} = [cos(angle)*ones(numel(obj.GeoPlotPointsQuiverRotor{patch}{1, 1}), 1), sin(angle)*ones(numel(obj.GeoPlotPointsQuiverRotor{patch}{2, 1}), 1)];
        end

        % Reset Phases
        function resetPhases(obj)
            obj.Phases=struct('Patches', {}, 'Type', {}, 'Slot', {}, 'Area', {}, 'RHS', {});
        end
        % Add a phase to 
        % type: "U" "V" or "W"
        % slot: -1 or 1
        function addPhase(obj, type, patches, slot)
            pos = numel(obj.Phases) + 1;
            obj.Phases(pos).Type = type;
            obj.Phases(pos).Patches = patches;
            obj.Phases(pos).Slot = slot;
            obj.Phases(pos).Area = sum(obj.PatchAreaStator(patches));
            obj.Phases(pos).RHS = op_f_v_mp_eval(obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MpMeshEvalSt, patches);
        end

        function setCurrentOptions(obj, applicationCurrent, polePairs, numWindings, phase0)
            obj.ApplicationCurrent = applicationCurrent;
            obj.PolePairs = polePairs;
            obj.NumWindings = numWindings;
            obj.Phase0 = deg2rad(phase0);            
        end

        function setCurrent(obj, iAapp, angle)
            angleRad = deg2rad(angle);
            obj.RHSStator = zeros(size(obj.RHSStator)); %reset old currents
            obj.dRHSStator = zeros(size(obj.RHSStator)); % derivatives wrt operating angle

            for iPhase = 1:numel(obj.Phases)
                switch obj.Phases(iPhase).Type
                    case "U"
                        Current = iAapp * sin (obj.Phase0 + obj.PolePairs * angleRad);
                        dCurrent = iAapp * cos (obj.Phase0 + obj.PolePairs * angleRad);
                    case "V"
                        Current = iAapp * sin (obj.Phase0 + obj.PolePairs * angleRad - 2/3*pi);
                        dCurrent = iAapp * cos (obj.Phase0 + obj.PolePairs * angleRad - 2/3*pi);
                    case "W"
                        Current = iAapp * sin (obj.Phase0 + obj.PolePairs * angleRad + 2/3*pi);
                        dCurrent = iAapp * cos (obj.Phase0 + obj.PolePairs * angleRad + 2/3*pi);
                    otherwise
                        error("Unknown phase type");
                end
                obj.RHSStator = obj.RHSStator ...
                    + obj.Phases(iPhase).RHS * Current * obj.NumWindings * obj.Phases(iPhase).Slot / obj.Phases(iPhase).Area;
                obj.dRHSStator = obj.dRHSStator ...
                    + obj.Phases(iPhase).RHS * dCurrent * obj.NumWindings * obj.Phases(iPhase).Slot / obj.Phases(iPhase).Area;
            end
        end
        
        %% Solve for the magnetic potential, linear case with antiperiodic BC
        function solveMagneticPotentialLinear(obj)
%             obj.generateStiffnessMatrices();
            KRotor = obj.StiffMatRotor;
            % set antiperiodicity
            KRotor(obj.APdofsRotorLeft, :) = KRotor(obj.APdofsRotorLeft, :) - KRotor(obj.APdofsRotorRight, :);
            KRotor(obj.APdofsRotorRight, :) = 0;
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorLeft)  = eye (numel (obj.APdofsRotorRight));
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorRight) = eye (numel (obj.APdofsRotorRight));

            KStator = obj.StiffMatStator;
            % set antiperiodicity
            KStator(obj.APdofsStatorLeft, :) = KStator(obj.APdofsStatorLeft, :) - KStator(obj.APdofsStatorRight, :);
            KStator(obj.APdofsStatorRight, :) = 0;
            KStator(obj.APdofsStatorRight, obj.APdofsStatorLeft)  = eye (numel (obj.APdofsStatorRight));
            KStator(obj.APdofsStatorRight, obj.APdofsStatorRight) = eye (numel (obj.APdofsStatorRight));

            % Assemble coupled rotor/stator system
            K = [KRotor(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor), sparse(numel(obj.IndependentDOFsRotor), numel(obj.IndependentDOFsStator));
                    sparse(numel(obj.IndependentDOFsStator), numel(obj.IndependentDOFsRotor)), KStator(obj.IndependentDOFsStator, obj.IndependentDOFsStator)];
            B = [-obj.CouplingMatrixRotor(obj.IndependentDOFsRotor, :) ; obj.CouplingMatrixStator(obj.IndependentDOFsStator, :)];

            A = [K, B;
                 B', sparse(obj.NumberDOFCoupling, obj.NumberDOFCoupling)];
            
            RHS = [obj.RHSRotor(obj.IndependentDOFsRotor); obj.RHSStator(obj.IndependentDOFsStator); zeros(obj.NumberDOFCoupling, 1)];
            % solve linear system
            u = A\RHS;

            obj.MagneticPotentialRotor = zeros(obj.NumberDOFRotor, 1);
            obj.MagneticPotentialStator = zeros(obj.NumberDOFStator, 1);
            % assign solution to rotor/stator/coupling
            obj.MagneticPotentialRotor(obj.IndependentDOFsRotor) = real(u(1:numel(obj.IndependentDOFsRotor)));
            obj.MagneticPotentialStator(obj.IndependentDOFsStator) = real(u(numel(obj.IndependentDOFsRotor)+1 : numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator)));
            obj.CouplingValues = real(u(numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator)+1 : end));

        end

        function solveMagneticPotentialNewton(obj)
            % If has not been solved yet, use linear solution as initial guess
            if all(obj.MagneticPotentialRotor==0)
                obj.solveMagneticPotentialLinear();
            end
            [urt, ust, ucouple, ~] = obj.newtonSolver(obj.MagneticPotentialRotor, obj.MagneticPotentialStator, ...
                obj.CouplingValues, obj.MultipatchSpaceRotor, obj.MpSpaceEvalRt, obj.MultipatchMeshRotor, obj.MpMeshEvalRt, ...
                obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MultipatchMeshStator, obj.MpMeshEvalSt);
            obj.MagneticPotentialRotor = urt;
            obj.MagneticPotentialStator = ust;
            obj.CouplingValues = ucouple;
        end

        function solveMagneticPotentialFixpoint(obj)
            % If has not been solved yet, use linear solution as initial guess
            if all(obj.MagneticPotentialRotor==0)
                obj.solveMagneticPotentialLinear();
            end
            [urt, ust, ucouple] = obj.fixpointSolver(obj.MagneticPotentialRotor, obj.MagneticPotentialStator, ...
                obj.CouplingValues, obj.MultipatchSpaceRotor, obj.MpSpaceEvalRt, obj.MultipatchMeshRotor, obj.MpMeshEvalRt, ...
                obj.MultipatchSpaceStator, obj.MpSpaceEvalSt, obj.MultipatchMeshStator, obj.MpMeshEvalSt);
            obj.MagneticPotentialRotor = urt;
            obj.MagneticPotentialStator = ust;
            obj.CouplingValues = ucouple;
        end

        function [urt, ust, ucouple, A] = linearSolver(obj, StiffMatRotor1, StiffMatStator1)
            KRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
            KStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
            for iMat = 1:numel(obj.Materials)
                KRotor = KRotor + obj.StiffMatsRotor{iMat};
                KStator = KStator + obj.StiffMatsStator{iMat};
            end
            % set antiperiodicity
            KRotor(obj.APdofsRotorLeft, :) = KRotor(obj.APdofsRotorLeft, :) - KRotor(obj.APdofsRotorRight, :);
            KRotor(obj.APdofsRotorRight, :) = 0;
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorLeft)  = eye (numel (obj.APdofsRotorRight));
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorRight) = eye (numel (obj.APdofsRotorRight));
            % set antiperiodicity
            KStator(obj.APdofsStatorLeft, :) = KStator(obj.APdofsStatorLeft, :) - KStator(obj.APdofsStatorRight, :);
            KStator(obj.APdofsStatorRight, :) = 0;
            KStator(obj.APdofsStatorRight, obj.APdofsStatorLeft)  = eye (numel (obj.APdofsStatorRight));
            KStator(obj.APdofsStatorRight, obj.APdofsStatorRight) = eye (numel (obj.APdofsStatorRight));

            % Assemble coupled rotor/stator system
            K = [KRotor(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor), sparse(numel(obj.IndependentDOFsRotor), numel(obj.IndependentDOFsStator));
                    sparse(numel(obj.IndependentDOFsStator), numel(obj.IndependentDOFsRotor)), KStator(obj.IndependentDOFsStator, obj.IndependentDOFsStator)];
            B = [-obj.CouplingMatrixRotor(obj.IndependentDOFsRotor, :) ; obj.CouplingMatrixStator(obj.IndependentDOFsStator, :)];

            A = [K, B;
                 B', sparse(obj.NumberDOFCoupling, obj.NumberDOFCoupling)];
            
            RHS = [obj.RHSRotor(obj.IndependentDOFsRotor); obj.RHSStator(obj.IndependentDOFsStator); zeros(obj.NumberDOFCoupling, 1)];

            % solve linear system
            u = A\RHS;

            urt = zeros(obj.NumberDOFRotor, 1);
            ust = zeros(obj.NumberDOFStator, 1);
            % assign solution to rotor/stator/coupling
            urt(obj.IndependentDOFsRotor) = real(u(1:numel(obj.IndependentDOFsRotor)));
            ust(obj.IndependentDOFsStator) = real(u(numel(obj.IndependentDOFsRotor)+1 : numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator)));
            ucouple = real(u(numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator)+1 : end));

        end

        function [urt, ust, ucouple, J] = newtonSolver(obj, urt, ust, ucouple, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, verbose)
            if nargin == 12
                verbose = true;
            end
            % Initial Evaluation
            u0 = [urt(obj.IndependentDOFsRotor); ust(obj.IndependentDOFsStator); ucouple];
            % update stiffness matrices
            KRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
            KStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
            for iMat = 1:numel(obj.Materials)
                if obj.Materials(iMat).Material.getType() == "nonlinear"
                    KRotor = KRotor + op_gradu_nu_gradv_mp_eval(sp_rt, sp_rt_eval, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, urt, obj.Materials(iMat).Material , obj.Materials(iMat).PatchesRotor);
                    KStator = KStator + op_gradu_nu_gradv_mp_eval(sp_st, sp_st_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, ust, obj.Materials(iMat).Material , obj.Materials(iMat).PatchesStator);
                else
                    KRotor = KRotor + obj.StiffMatsRotor{iMat};
                    KStator = KStator + obj.StiffMatsStator{iMat};
                end
            end
            % Antiperiodic boundary conditions
            % Rotor
            KRotor(obj.APdofsRotorLeft, :) = KRotor(obj.APdofsRotorLeft, :) - KRotor(obj.APdofsRotorRight, :);
            KRotor(obj.APdofsRotorRight, :) = 0;
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorLeft)  = eye (numel (obj.APdofsRotorRight));
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorRight) = eye (numel (obj.APdofsRotorRight));
            % Stator
            KStator(obj.APdofsStatorLeft, :) = KStator(obj.APdofsStatorLeft, :) - KStator(obj.APdofsStatorRight, :);
            KStator(obj.APdofsStatorRight, :) = 0;
            KStator(obj.APdofsStatorRight, obj.APdofsStatorLeft)  = eye (numel (obj.APdofsStatorRight));
            KStator(obj.APdofsStatorRight, obj.APdofsStatorRight) = eye (numel (obj.APdofsStatorRight));

            % Assemble coupled rotor/stator system
            K = [KRotor(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor), sparse(numel(obj.IndependentDOFsRotor), numel(obj.IndependentDOFsStator));
                    sparse(numel(obj.IndependentDOFsStator), numel(obj.IndependentDOFsRotor)), KStator(obj.IndependentDOFsStator, obj.IndependentDOFsStator)];
            B = [-obj.CouplingMatrixRotor(obj.IndependentDOFsRotor, :) ; obj.CouplingMatrixStator(obj.IndependentDOFsStator, :)];

            A = [K, B;
                 B', sparse(obj.NumberDOFCoupling, obj.NumberDOFCoupling)];
            f = [obj.RHSRotor(obj.IndependentDOFsRotor); obj.RHSStator(obj.IndependentDOFsStator); zeros(numel(ucouple), 1)];

            r0 = A*u0-f;
            residual0 = r0'*r0;
            
            % Newton Iteration
            for NewtonIter = 1:50
                % Calculate Jacobian Matrix
                K_rt_prime = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
                K_st_prime = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
                for iMat = 1:numel(obj.Materials)
                    if obj.Materials(iMat).Material.getType() == "nonlinear"
                        K_rt_prime = K_rt_prime + op_dK_du_times_u_mp_eval(sp_rt, sp_rt_eval, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, urt, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesRotor);
                        K_st_prime = K_st_prime + op_dK_du_times_u_mp_eval(sp_st, sp_st_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, ust, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesStator);
                    end
                end
                % K_rt_prime = op_dK_du_times_u_mp_eval(sp_rt, sp_rt_eval, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, urt, obj.Material, obj.PatchesRotor.Iron);
                % Stator
                % K_st_prime = op_dK_du_times_u_mp_eval(sp_st, sp_st_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, ust, obj.Material, obj.PatchesStator.Iron);
                % Antiperiodic boundary conditions
                K_rt_prime(obj.APdofsRotorLeft', :) = K_rt_prime(obj.APdofsRotorLeft', :) - K_rt_prime(obj.APdofsRotorRight', :);
                K_rt_prime(obj.APdofsRotorRight', :) = 0;
                K_st_prime(obj.APdofsStatorLeft', :) = K_st_prime(obj.APdofsStatorLeft', :) - K_st_prime(obj.APdofsStatorRight', :);
                K_st_prime(obj.APdofsStatorRight', :) = 0;
                % Dirichlet boundary conditions
                K_rt_prime_shrink = K_rt_prime(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor);
                K_st_prime_shrink = K_st_prime(obj.IndependentDOFsStator, obj.IndependentDOFsStator);
                % Jacobian matrix
                J = A + [K_rt_prime_shrink, sparse(size(K_rt_prime_shrink, 2), size(K_st_prime_shrink, 1)+numel(ucouple));
                         sparse(size(K_st_prime_shrink, 1), size(K_rt_prime_shrink, 2)), K_st_prime_shrink, sparse(size(K_st_prime_shrink, 1), numel(ucouple));
                         sparse(numel(ucouple), size(A, 2))];

                w = J\r0;

                tau = 1;
                % Linesearch
                Lmax = 10;
                for L = 1:Lmax
                    u1 = u0 - tau*w;
                    urt(obj.IndependentDOFsRotor) = u1(1:numel(obj.IndependentDOFsRotor));
                    ust(obj.IndependentDOFsStator) = u1(numel(obj.IndependentDOFsRotor)+(1:numel(obj.IndependentDOFsStator)));
                    ucouple = u1(numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator)+1:end);

                    % update stiffness matrices
                    KRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
                    KStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
                    for iMat = 1:numel(obj.Materials)
                        if obj.Materials(iMat).Material.getType() == "nonlinear"
                            KRotor = KRotor + op_gradu_nu_gradv_mp_eval(sp_rt, sp_rt_eval, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, urt, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesRotor);
                            KStator = KStator + op_gradu_nu_gradv_mp_eval(sp_st, sp_st_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, ust, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesStator);
                        else
                            KRotor = KRotor + obj.StiffMatsRotor{iMat};
                            KStator = KStator + obj.StiffMatsStator{iMat};
                        end
                    end
                    % Antiperiodic boundary conditions
                    % Rotor
                    KRotor(obj.APdofsRotorLeft, :) = KRotor(obj.APdofsRotorLeft, :) - KRotor(obj.APdofsRotorRight, :);
                    KRotor(obj.APdofsRotorRight, :) = 0;
                    KRotor(obj.APdofsRotorRight, obj.APdofsRotorLeft)  = eye (numel (obj.APdofsRotorRight));
                    KRotor(obj.APdofsRotorRight, obj.APdofsRotorRight) = eye (numel (obj.APdofsRotorRight));
                    % Stator
                    KStator(obj.APdofsStatorLeft, :) = KStator(obj.APdofsStatorLeft, :) - KStator(obj.APdofsStatorRight, :);
                    KStator(obj.APdofsStatorRight, :) = 0;
                    KStator(obj.APdofsStatorRight, obj.APdofsStatorLeft)  = eye (numel (obj.APdofsStatorRight));
                    KStator(obj.APdofsStatorRight, obj.APdofsStatorRight) = eye (numel (obj.APdofsStatorRight));

                    % Assemble coupled rotor/stator system
                    K = [KRotor(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor), sparse(numel(obj.IndependentDOFsRotor), numel(obj.IndependentDOFsStator));
                            sparse(numel(obj.IndependentDOFsStator), numel(obj.IndependentDOFsRotor)), KStator(obj.IndependentDOFsStator, obj.IndependentDOFsStator)];
                    B = [-obj.CouplingMatrixRotor(obj.IndependentDOFsRotor, :) ; obj.CouplingMatrixStator(obj.IndependentDOFsStator, :)];
        
                    A = [K, B;
                         B', sparse(obj.NumberDOFCoupling, obj.NumberDOFCoupling)];
                    f = [obj.RHSRotor(obj.IndependentDOFsRotor); obj.RHSStator(obj.IndependentDOFsStator); zeros(numel(obj.CouplingValues), 1)];

                    r1 = (A*u1 - f);
                    residual1 = r1'*r1;

                    if (residual1 < residual0)
                        if verbose == true
                            disp(['NewtonIter: ' 9 num2str(NewtonIter) 9 'Residual: ' 9 num2str(residual1) 9 'Linesearch ' num2str(L)])
                        end
                        u0 = u1;
                        r0 = r1;
                        residual0 = residual1;
                        break;
                    else
                        tau = tau/2;
                    end
                    if L == Lmax
                        % switch to fixpoint for some iterations
                        disp("Max Iteration in Newton Linesearch reached, switched to fixpoint!")
                        [urt, ust, ucouple, r0] = obj.fixpointSolver(urt, ust, ucouple, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, verbose);
                        u0 = [urt(obj.IndependentDOFsRotor); ust(obj.IndependentDOFsStator); ucouple];
                        residual0 = r0'*r0 ;
                        if residual0 < obj.NewtonResidual
                            return
                        end
                    end
                end
                if residual1 < obj.NewtonResidual
                    break;
                end
            end
        end

        function [urt, ust, ucouple, r0] = fixpointSolver(obj, urt, ust, ucouple, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, verbose)
            if nargin == 12
                verbose = true;
            end
            % Initial Evaluation
            u0 = [urt(obj.IndependentDOFsRotor); ust(obj.IndependentDOFsStator); ucouple];
            % update stiffness matrices
            KRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
            KStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
            for iMat = 1:numel(obj.Materials)
                if obj.Materials(iMat).Material.getType() == "nonlinear"
                    KRotor = KRotor + op_gradu_nu_gradv_mp_eval(sp_rt, sp_rt_eval, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, urt, obj.Materials(iMat).Material , obj.Materials(iMat).PatchesRotor);
                    KStator = KStator + op_gradu_nu_gradv_mp_eval(sp_st, sp_st_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, ust, obj.Materials(iMat).Material , obj.Materials(iMat).PatchesStator);
                else
                    KRotor = KRotor + obj.StiffMatsRotor{iMat};
                    KStator = KStator + obj.StiffMatsStator{iMat};
                end
            end
            % Antiperiodic boundary conditions
            % Rotor
            KRotor(obj.APdofsRotorLeft, :) = KRotor(obj.APdofsRotorLeft, :) - KRotor(obj.APdofsRotorRight, :);
            KRotor(obj.APdofsRotorRight, :) = 0;
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorLeft)  = eye (numel (obj.APdofsRotorRight));
            KRotor(obj.APdofsRotorRight, obj.APdofsRotorRight) = eye (numel (obj.APdofsRotorRight));
            % Stator
            KStator(obj.APdofsStatorLeft, :) = KStator(obj.APdofsStatorLeft, :) - KStator(obj.APdofsStatorRight, :);
            KStator(obj.APdofsStatorRight, :) = 0;
            KStator(obj.APdofsStatorRight, obj.APdofsStatorLeft)  = eye (numel (obj.APdofsStatorRight));
            KStator(obj.APdofsStatorRight, obj.APdofsStatorRight) = eye (numel (obj.APdofsStatorRight));

            % Assemble coupled rotor/stator system
            K = [KRotor(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor), sparse(numel(obj.IndependentDOFsRotor), numel(obj.IndependentDOFsStator));
                    sparse(numel(obj.IndependentDOFsStator), numel(obj.IndependentDOFsRotor)), KStator(obj.IndependentDOFsStator, obj.IndependentDOFsStator)];
            B = [-obj.CouplingMatrixRotor(obj.IndependentDOFsRotor, :) ; obj.CouplingMatrixStator(obj.IndependentDOFsStator, :)];

            A = [K, B;
                 B', sparse(obj.NumberDOFCoupling, obj.NumberDOFCoupling)];
            f = [obj.RHSRotor(obj.IndependentDOFsRotor); obj.RHSStator(obj.IndependentDOFsStator); zeros(numel(obj.CouplingValues), 1)];

            r0 = A*u0-f;
            residual0 = r0'*r0;
            
            % Newton Iteration
            for FixpointIter = 1:10
                w = u0 - A\f;

                tau = 1/4;
                % Linesearch
                for L = 1:10
                    u1 = u0 - tau*w;
                    urt(obj.IndependentDOFsRotor) = u1(1:numel(obj.IndependentDOFsRotor));
                    ust(obj.IndependentDOFsStator) = u1(numel(obj.IndependentDOFsRotor)+(1:numel(obj.IndependentDOFsStator)));
                    ucouple = u1(numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator)+1:end);

                    % update stiffness matrices
                    KRotor = sparse(obj.NumberDOFRotor, obj.NumberDOFRotor);
                    KStator = sparse(obj.NumberDOFStator, obj.NumberDOFStator);
                    for iMat = 1:numel(obj.Materials)
                        if obj.Materials(iMat).Material.getType() == "nonlinear"
                            KRotor = KRotor + op_gradu_nu_gradv_mp_eval(sp_rt, sp_rt_eval, sp_rt, sp_rt_eval, msh_rt, msh_rt_eval, urt, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesRotor);
                            KStator = KStator + op_gradu_nu_gradv_mp_eval(sp_st, sp_st_eval, sp_st, sp_st_eval, msh_st, msh_st_eval, ust, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesStator);
                        else
                            KRotor = KRotor + obj.StiffMatsRotor{iMat};
                            KStator = KStator + obj.StiffMatsStator{iMat};
                        end
                    end
                    % Antiperiodic boundary conditions
                    % Rotor
                    KRotor(obj.APdofsRotorLeft, :) = KRotor(obj.APdofsRotorLeft, :) - KRotor(obj.APdofsRotorRight, :);
                    KRotor(obj.APdofsRotorRight, :) = 0;
                    KRotor(obj.APdofsRotorRight, obj.APdofsRotorLeft)  = eye (numel (obj.APdofsRotorRight));
                    KRotor(obj.APdofsRotorRight, obj.APdofsRotorRight) = eye (numel (obj.APdofsRotorRight));
                    % Stator
                    KStator(obj.APdofsStatorLeft, :) = KStator(obj.APdofsStatorLeft, :) - KStator(obj.APdofsStatorRight, :);
                    KStator(obj.APdofsStatorRight, :) = 0;
                    KStator(obj.APdofsStatorRight, obj.APdofsStatorLeft)  = eye (numel (obj.APdofsStatorRight));
                    KStator(obj.APdofsStatorRight, obj.APdofsStatorRight) = eye (numel (obj.APdofsStatorRight));

                    % Assemble coupled rotor/stator system
                    K = [KRotor(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor), sparse(numel(obj.IndependentDOFsRotor), numel(obj.IndependentDOFsStator));
                            sparse(numel(obj.IndependentDOFsStator), numel(obj.IndependentDOFsRotor)), KStator(obj.IndependentDOFsStator, obj.IndependentDOFsStator)];
                    B = [-obj.CouplingMatrixRotor(obj.IndependentDOFsRotor, :) ; obj.CouplingMatrixStator(obj.IndependentDOFsStator, :)];
        
                    A = [K, B;
                         B', sparse(obj.NumberDOFCoupling, obj.NumberDOFCoupling)];
                    f = [obj.RHSRotor(obj.IndependentDOFsRotor); obj.RHSStator(obj.IndependentDOFsStator); zeros(numel(obj.CouplingValues), 1)];

                    r1 = (A*u1 - f);
                    residual1 = r1'*r1;

                    
                    if (residual1 < residual0)
                        if verbose == true
                            disp(['FixpointIter: ' 9 num2str(FixpointIter) 9 'Residual: ' 9 num2str(residual1) 9 'Linesearch ' num2str(L)])
                        end
                        u0 = u1;
                        r0 = r1;
                        residual0 = residual1;      
%                         obj.plotBResulting();
                        break;
                    else
                        tau = tau/2;
                    end
                end
                if residual1 < obj.NewtonResidual
                    break;
                end
%                 obj.plotBResulting();
            end

        end

        %% Torque calculation function
        function T = calcTorqueCoupling(obj)
            T = - obj.MagneticPotentialStator' * obj.CouplingMatrixStatorInit * obj.RotationMatrixDer * obj.CouplingValues * obj.MachineProperties.Length;
        end

        function T = calcTorqueBrBtRotor(obj)
            T = 0;
            for bnd = 1:numel(obj.BoundariesCouplingRotor)
                patch = obj.BoundariesRotor(obj.BoundariesCouplingRotor(bnd)).patches;
                u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{patch});
                T = T +  1/obj.MuVac * obj.MachineProperties.Length * ...
                op_Br_Bt_dGamma(obj.CouplingSpacesRotorEval{bnd}, obj.CouplingMeshesRotor{bnd}, u_loc);
            end
        end

        function T = calcTorqueBrBtStator(obj)
            T = 0;
            for bnd = 1:numel(obj.BoundariesCouplingStator)
                patch = obj.BoundariesStator(obj.BoundariesCouplingStator(bnd)).patches;
                u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{patch});
                T = T +  1/obj.MuVac * obj.MachineProperties.Length * ...
                op_Br_Bt_dGamma(obj.CouplingSpacesStatorEval{bnd}, obj.CouplingMeshesStator{bnd}, u_loc);
            end
        end

        function calcMagneticForces(obj)
            obj.ForceDensityRotor = zeros(obj.MultipatchSpaceRotor.ndof, 2);
            obj.ForceDensityStator = zeros(obj.MultipatchSpaceStator.ndof, 2);
            for iMat = 1:numel(obj.Materials)
                obj.ForceDensityRotor = obj.ForceDensityRotor + op_maxwell_force_mp(obj.MultipatchSpaceRotor, obj.MultipatchMeshRotor, obj.MagneticPotentialRotor, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesRotor);
                obj.ForceDensityStator = obj.ForceDensityStator + op_maxwell_force_mp(obj.MultipatchSpaceStator, obj.MultipatchMeshStator, obj.MagneticPotentialStator, obj.Materials(iMat).Material, obj.Materials(iMat).PatchesStator);
            end
        end

        function plotMagneticForcesEXPERIMENTAL(obj)
            % Quivers
                Fx_vals = [];
                Fy_vals = [];
                X_vals = [];
                Y_vals =[];
                for iPatch = 1:obj.NumberPatchesRotor
                    fx_loc = obj.ForceDensityRotor(obj.MultipatchSpaceRotor.gnum{iPatch}, 1);
                    fy_loc = obj.ForceDensityRotor(obj.MultipatchSpaceRotor.gnum{iPatch}, 2);
                    Fx = sum(obj.RotorShapeValuesQuiver{iPatch}{1}.* reshape(fx_loc, 1, 1, []), 3);
                    Fy = sum(obj.RotorShapeValuesQuiver{iPatch}{1}.* reshape(fy_loc, 1, 1, []), 3);
                    %
                    x = obj.GeoPlotPointsQuiverRotor{iPatch}{1, 1};
                    y = obj.GeoPlotPointsQuiverRotor{iPatch}{2, 1};

                    X_vals = cat(2, X_vals, reshape(x, 1, []));
                    Y_vals = cat(2, Y_vals, reshape(y, 1, []));
                    Fx_vals = cat(2, Fx_vals, reshape(Fx, 1, []));
                    Fy_vals = cat(2, Fy_vals, reshape(Fy, 1, []));
                end
                % Rotate Rotor quivers
                Fxtemp = Fx_vals*cos(obj.RotationAngle) - Fy_vals*sin(obj.RotationAngle);
                Fytemp = Fy_vals*cos(obj.RotationAngle) + Fx_vals*sin(obj.RotationAngle);
                Fx_vals = Fxtemp;
                Fy_vals = Fytemp;

                for iPatch = 1:obj.NumberPatchesStator
                    fx_loc = obj.ForceDensityStator(obj.MultipatchSpaceStator.gnum{iPatch}, 1);
                    fy_loc = obj.ForceDensityStator(obj.MultipatchSpaceStator.gnum{iPatch}, 2);
                    Fx = sum(obj.StatorShapeValuesQuiver{iPatch}{1}.* reshape(fx_loc, 1, 1, []), 3);
                    Fy = sum(obj.StatorShapeValuesQuiver{iPatch}{1}.* reshape(fy_loc, 1, 1, []), 3);
                    %
                    x = obj.GeoPlotPointsQuiverStator{iPatch}{1, 1};
                    y = obj.GeoPlotPointsQuiverStator{iPatch}{2, 1};

                    X_vals = cat(2, X_vals, reshape(x, 1, []));
                    Y_vals = cat(2, Y_vals, reshape(y, 1, []));
                    Fx_vals = cat(2, Fx_vals, reshape(Fx, 1, []));
                    Fy_vals = cat(2, Fy_vals, reshape(Fy, 1, []));
                end

                quiver(X_vals, Y_vals, Fx_vals, Fy_vals, 'color', 'k');
%                 clims = [min(cols), max(cols)];
%                 caxis(clims);
                c = colorbar();
                c.TickLabelInterpreter = 'latex';
                c.Label.Interpreter = 'latex';
                c.Label.String = 'Magnetic flux density (T)';
                c.Label.S
                view([0, 90]);
                axis equal;
                xlim(obj.XLimits);
                ylim(obj.YLimits);
        end


        %% Calculate the Shape Values and Gradients at Plot Points in Advance
        function calcShapeValues(obj)
            disp("Calculating Shape Values and Derivatives at Plot Points...")
            for iPatch = 1:obj.NumberPatchesRotor
                obj.RotorShapeValues{iPatch} = shp_eval1(obj.MultipatchSpaceRotor.sp_patch{iPatch}, obj.GeometryRotor(iPatch), obj.ParamPlotPointsRotor{iPatch}, {'value', 'gradient'});
                obj.RotorShapeValuesQuiver{iPatch} = shp_eval1(obj.MultipatchSpaceRotor.sp_patch{iPatch}, obj.GeometryRotor(iPatch), obj.ParamPlotPointsQuiverRotor{iPatch}, {'value', 'gradient'});
            end
            for iPatch = 1:obj.NumberPatchesStator
                obj.StatorShapeValues{iPatch} = shp_eval1(obj.MultipatchSpaceStator.sp_patch{iPatch}, obj.GeometryStator(iPatch), obj.ParamPlotPointsStator{iPatch}, {'value', 'gradient'});
                obj.StatorShapeValuesQuiver{iPatch} = shp_eval1(obj.MultipatchSpaceStator.sp_patch{iPatch}, obj.GeometryStator(iPatch), obj.ParamPlotPointsQuiverStator{iPatch}, {'value', 'gradient'});
            end
        end

        % Plots -------------------------
        %% Plot the geometry
        function plotGeometry(obj)
            obj.calculateGeoPlotPoints();
            obj.FigureGeometry = figure('Name', 'Geometry Plot', 'NumberTitle', 'off', 'Position', [0 0 800 500]);
            hold on
            axis equal
            obj.plotPatches();
%             obj.plotPatchBoundaries();
            obj.plotMaterialBoundaries();
            set(obj.FigureGeometry.CurrentAxes, "XLim", obj.XLimits)
            set(obj.FigureGeometry.CurrentAxes, "YLim", obj.YLimits)
            obj.plotRemanenceQuiver();
        end

        function plotPatches(obj)
            for iMat = 1:numel(obj.Materials)
                for patch = obj.Materials(iMat).PatchesRotor
                    obj.FigureGeometryRotor{patch}=surf(obj.GeoPlotPointsRotor{patch}{1, 1}, obj.GeoPlotPointsRotor{patch}{2, 1}, zeros(size(obj.GeoPlotPointsRotor{patch}{1, 1})), ...
                        "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.PlotFaceAlpha);
                end
                for patch = obj.Materials(iMat).PatchesStator
                    obj.FigureGeometryStator{patch}=surf(obj.GeoPlotPointsStator{patch}{1, 1}, obj.GeoPlotPointsStator{patch}{2, 1}, zeros(size(obj.GeoPlotPointsStator{patch}{1, 1})), ...
                        "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.PlotFaceAlpha);
                end
            end
        end

        function plotPatchBoundaries(obj)
            % Rotor
            for iBnd = 1:numel(obj.BoundariesRotor)
                patch = obj.BoundariesRotor(iBnd).patches;
                side = obj.BoundariesRotor(iBnd).faces;
                obj.FigureBoundariesRotor{iBnd} = plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
            end
            for iInt = 1:numel(obj.InterfacesRotor)
                patch = obj.InterfacesRotor(iInt).patch1;
                side = obj.InterfacesRotor(iInt).side1;
                obj.FigureInterfacesRotor{iInt} = plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
            end
            % Stator
            for iBnd = 1:numel(obj.BoundariesStator)
                patch = obj.BoundariesStator(iBnd).patches;
                side = obj.BoundariesStator(iBnd).faces;
                obj.FigureBoundariesStator{iBnd} = plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
            end
            for iInt = 1:numel(obj.InterfacesStator)
                patch = obj.InterfacesStator(iInt).patch1;
                side = obj.InterfacesStator(iInt).side1;
                obj.FigureInterfacesStator{iInt} = plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
            end
        end

        function plotMaterialBoundaries(obj)
            % Rotor
            % Interfaces between different Materials
            for iInt = 1:numel(obj.InterfacesRotor)
                patch1 = obj.InterfacesRotor(iInt).patch1;
                side1 = obj.InterfacesRotor(iInt).side1;
                patch2 = obj.InterfacesRotor(iInt).patch2;
                % Only plot boundaries for interfaces with different materials
                plotInterface = true;
                for iMat = 1:numel(obj.Materials)
                    if isempty(obj.Materials(iMat).PatchesRotor)
                        continue
                    end
                    if any(obj.Materials(iMat).PatchesRotor == patch1) && any(obj.Materials(iMat).PatchesRotor == patch2)
                        plotInterface = false;
                    end
                end
                % Check magnets and set true again if different remanence direction
                % TB REVISED!
                for iRem = 1:numel(obj.Remanence)
                    if any (obj.Remanence(iRem).Patch == patch1)
                        for jRem = 1:numel(obj.Remanence)
                            if any (obj.Remanence(jRem).Patch == patch2)
                                if iRem == jRem
                                    continue
                                end
                                if obj.Remanence(iRem).InitAngle ~= obj.Remanence(jRem).InitAngle
                                    plotInterface = true;
                                end
                            end
                        end 
                    end
                end 
                % Plot boundaries
                if plotInterface == true
                    plot(obj.GeoPlotPointsRotor{patch1}{1, side1+1}, obj.GeoPlotPointsRotor{patch1}{2, side1+1}, "Color", obj.PlotColor.Contour, "LineWidth", 1.5);
                end
            end
            % Outside boundary
            for iInt = 1:numel(obj.BoundariesRotor)
                patch = obj.BoundariesRotor(iInt).patches;
                side = obj.BoundariesRotor(iInt).faces;
                plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Contour, "LineWidth", 1.5);
            end

            % Stator
            % Interfaces between different Materials
            for iInt = 1:numel(obj.InterfacesStator)
                patch1 = obj.InterfacesStator(iInt).patch1;
                side1 = obj.InterfacesStator(iInt).side1;
                patch2 = obj.InterfacesStator(iInt).patch2;
                % Only plot boundaries for interfaces with different materials
                plotInterface = true;
                for iMat = 1:numel(obj.Materials)
                    if isempty(obj.Materials(iMat).PatchesStator)
                        continue
                    end
                    if any(obj.Materials(iMat).PatchesStator == patch1) && any(obj.Materials(iMat).PatchesStator == patch2)
                        plotInterface = false;
                    end
                end
                % Check if interface is between two phases, then set plot true again
                for iPhase = 1:numel(obj.Phases)
                    if any (obj.Phases(iPhase).Patches == patch1)
                        for jPhase = 1:numel(obj.Phases)
                            if iPhase == jPhase
                                continue
                            end
                            if any (obj.Phases(jPhase).Patches == patch2)
                                plotInterface = true;
                            end
                        end 
                    end
                end 
                % Plot boundaries
                if plotInterface == true
                    plot(obj.GeoPlotPointsStator{patch1}{1, side1+1}, obj.GeoPlotPointsStator{patch1}{2, side1+1}, "Color", obj.PlotColor.Contour, "LineWidth", 1.5);
                end
            end
            % Outside boundary
            for iInt = 1:numel(obj.BoundariesStator)
                patch = obj.BoundariesStator(iInt).patches;
                side = obj.BoundariesStator(iInt).faces;
                plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Contour, "LineWidth", 1.5);
            end
        end

        % Plot the boundary conditions (Coupling/Dirichlet/Antiperiodic)
        function plotBoundariesRB(obj)
            % Rotor
            lw = 2;
            for iBnd = 1:numel(obj.BoundariesRotor)
                patch1 = obj.BoundariesRotor(iBnd).patches;
                side1 = obj.BoundariesRotor(iBnd).faces;
                if any(iBnd == obj.BoundariesDirichletRotor)
                    plot(obj.GeoPlotPointsRotor{patch1}{1, side1+1}, obj.GeoPlotPointsRotor{patch1}{2, side1+1}, "Color", TUDa_getColor('2b'), "LineWidth", lw);
                elseif any(iBnd == obj.BoundariesCouplingRotor)
                    plot(obj.GeoPlotPointsRotor{patch1}{1, side1+1}, obj.GeoPlotPointsRotor{patch1}{2, side1+1}, "Color", TUDa_getColor('9b'), "LineWidth", lw);
                elseif any(iBnd == obj.BoundariesAntiperiodicRotor)
                    plot(obj.GeoPlotPointsRotor{patch1}{1, side1+1}, obj.GeoPlotPointsRotor{patch1}{2, side1+1}, "Color", TUDa_getColor('5a'), "LineWidth", lw);
                else
                    warning("Not all boundaries are defined!")
                end
            end
            % Stator
            for iBnd = 1:numel(obj.BoundariesStator)
                patch1 = obj.BoundariesStator(iBnd).patches;
                side1 = obj.BoundariesStator(iBnd).faces;
                if any(iBnd == obj.BoundariesDirichletStator)
                    plot(obj.GeoPlotPointsStator{patch1}{1, side1+1}, obj.GeoPlotPointsStator{patch1}{2, side1+1}, "Color", TUDa_getColor('2b'), "LineWidth", lw);
                elseif any(iBnd == obj.BoundariesCouplingStator)
                    plot(obj.GeoPlotPointsStator{patch1}{1, side1+1}, obj.GeoPlotPointsStator{patch1}{2, side1+1}, "Color", TUDa_getColor('9b'), "LineWidth", lw);
                elseif any(iBnd == obj.BoundariesAntiperiodicStator)
                    plot(obj.GeoPlotPointsStator{patch1}{1, side1+1}, obj.GeoPlotPointsStator{patch1}{2, side1+1}, "Color", TUDa_getColor('5a'), "LineWidth", lw);
                else
                    warning("Not all boundaries are defined!")
                end
            end
        end
        
        % Plot the winding currents
        function plotWindingCurrent(obj)
            function plotInCurrent(x, y, r)
                p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
                p1 = nsidedpoly(1000, 'Center', [x y], 'Radius', r/5);
                plot(p, 'FaceColor', 'white', 'FaceAlpha', 1, 'LineWidth', 3)
                plot(p1, 'FaceColor', 'black', 'FaceAlpha', 1)
            end
            function plotOutCurrent(x, y, r, rot)
                p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
                pts = [-r, -r, r, r; -r, r -r, r]/(2^0.5);
                R = [cos(rot), -sin(rot); sin(rot), cos(rot)];
                pts = R*pts;
                plot(p, 'FaceColor', 'white', 'FaceAlpha', 1, 'LineWidth', 3)
                plot(x + [pts(1, 1), pts(1, 4) ], y + [pts(2, 1), pts(2, 4) ], 'LineWidth', 3, 'Color', 'black')
                plot(x + [pts(1, 2), pts(1, 3) ], y + [pts(2, 2), pts(2, 3) ], 'LineWidth', 3, 'Color', 'black')
            end
            for iPhase = 1:numel(obj.Phases)
                AreaMoment = 0;
                for j = 1:numel(obj.Phases(iPhase).Patches)
                    patch = obj.Phases(iPhase).Patches(j);
                    AreaMoment = AreaMoment + obj.PatchAreaStator(patch)*nrbeval(obj.GeometryStator(patch).nurbs, {0.5, 0.5});
                end
                Xmean = AreaMoment/obj.Phases(iPhase).Area;
                angle = atan2(Xmean(2), Xmean(1));
                if obj.Phases(iPhase).Slot == 1
                    plotInCurrent(Xmean(1), Xmean(2), obj.Phases(iPhase).Area^0.5/5);
                else
                    plotOutCurrent(Xmean(1), Xmean(2), obj.Phases(iPhase).Area^0.5/5, angle);
                end
            end
        end

        function plotControlPoints(obj, rotstat, drawText, indRot, indStat, markSize, col)
            if nargin <= 1
                rotstat = [1, 2];
            end
            if nargin <= 2
                drawText = false;
            end
            if nargin <= 3 || isempty(indRot)
                indRot = 1:size(obj.ControlPointsRotor, 1);
            end
            if nargin <= 4 || isempty(indStat)
                indStat = 1:size(obj.ControlPointsStator, 1);
            end
            if nargin <= 5 || isempty(markSize)
                markSize = 15;
            end
            if nargin <= 6
                col = TUDa_getColor("9b");
            end

            % Rotor
            if any(rotstat ==1)
                scatter(obj.ControlPointsRotor(indRot, 1), obj.ControlPointsRotor(indRot, 2), "filled", 'MarkerFaceColor', col, 'SizeData', markSize, Marker='o')
                if drawText
                text(obj.ControlPointsRotor(indRot, 1), obj.ControlPointsRotor(indRot, 2), string(indRot));
                end
            end
            % Stator
            if any(rotstat ==2)
                scatter(obj.ControlPointsStator(indStat, 1), obj.ControlPointsStator(indStat, 2), "filled", 'MarkerFaceColor', col, 'SizeData', markSize, Marker='o')
                if drawText
                text(obj.ControlPointsStator(indStat, 1), obj.ControlPointsStator(indStat, 2), string(indStat));
                end
            end
            
        end
        %% For testing, get boundary numbers
        function plotBoundaryNr(obj, rotStat)
            if nargin == 1 
                rotStat = [1, 2];
            end
            % Rotor
            if any(rotStat == 1)
                for Num = 1:numel(obj.BoundariesRotor)
                    x = nrbeval(obj.GeometryRotor(obj.BoundariesRotor(Num).patches).boundary(obj.BoundariesRotor(Num).faces).nurbs, 0.5);
                    text(x(1), x(2), num2str(Num), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
                end
            end
            % Stator
            if any(rotStat == 2)
                for Num = 1:numel(obj.BoundariesStator)
                    x = nrbeval(obj.GeometryStator(obj.BoundariesStator(Num).patches).boundary(obj.BoundariesStator(Num).faces).nurbs, 0.5);
                    text(x(1), x(2), num2str(Num), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
                end
            end
        end

        %% Plot the patch numbers
        % input: RotStat - 1: plot Rotor patch numbers; 
        %                  2: plot Stator patch numbers
        function plotPatchNr(obj, rotStat)
            if nargin == 1 
                rotStat = [1, 2];
            end
            % Rotor
            if any(rotStat == 1)
                for Num = 1:obj.NumberPatchesRotor
                    x = obj.MultipatchMeshRotor.msh_patch{Num}.map([0.5, 0.5]);
                    text(x(1), x(2), num2str(Num), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
                end
            end
            % Stator
            if any(rotStat == 2)
                for Num = 1:obj.NumberPatchesStator
                    x = obj.MultipatchMeshStator.msh_patch{Num}.map([0.5, 0.5]);
                    text(x(1), x(2), num2str(Num), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
                end
            end
        end
        
        %% Calculate all patch areas and Inertias
        function calculateAreaAndInertia(obj)
            disp("Calculating Patch Areas...");

            for iPatch = 1:obj.NumberPatchesRotor
                msh = msh_precompute(obj.MultipatchMeshRotor.msh_patch{iPatch});
                obj.PatchAreaRotor(iPatch) = sum(msh.jacdet .* msh.quad_weights, 'all');
                obj.PatchInertiaRotor(iPatch) = sum(msh.jacdet .* msh.quad_weights .* reshape((msh.geo_map(1, :, :).^2 + msh.geo_map(2, :, :).^2) , msh.nqn, msh.nel), 'all');

            end
            for iPatch = 1:obj.NumberPatchesStator
                msh = msh_precompute(obj.MultipatchMeshStator.msh_patch{iPatch});
                obj.PatchAreaStator(iPatch) = sum(msh.jacdet .* msh.quad_weights, 'all');
                obj.PatchInertiaStator(iPatch) = sum(msh.jacdet .* msh.quad_weights .* reshape((msh.geo_map(1, :, :).^2 + msh.geo_map(2, :, :).^2) , msh.nqn, msh.nel), 'all');
            end
        end

        %% Calculate the plotting points for the quivers, depending on size of each patch
        function calculateQuiverPlotPoints(obj)
            disp("Calculating Quiver Plot Points...")
            % Rotor
            A_ges_rotor = sum(obj.PatchAreaRotor);
            for iPatch = 1:obj.NumberPatchesRotor
                % lengths of sides L1: left+right, L2: top+bottom
                L1 = sum(msh_eval_boundary_side(obj.MultipatchMeshRotor.msh_patch{iPatch}, 1).element_size) + ...
                   sum(msh_eval_boundary_side(obj.MultipatchMeshRotor.msh_patch{iPatch}, 2).element_size); 
                L2 = sum(msh_eval_boundary_side(obj.MultipatchMeshRotor.msh_patch{iPatch}, 3).element_size) + ...
                   sum(msh_eval_boundary_side(obj.MultipatchMeshRotor.msh_patch{iPatch}, 4).element_size); 
                num_quivers_i = ceil(obj.PatchAreaRotor(iPatch)/A_ges_rotor*obj.NumberQuiversRotor);
                % Number quivers in x an y according to lengths of sides
                q1 = ceil((num_quivers_i* L2/L1)^0.5);
                q2 = ceil((num_quivers_i* L1/L2)^0.5);
                % Points with also half distance to boundary
                obj.ParamPlotPointsQuiverRotor{iPatch} = {linspace(0, 1-1/q1, q1)+0.5/q1, linspace(0, 1-1/q2, q2)+0.5/q2};
                % Points also directly on boundary
%                 obj.ParamPlotPointsQuiverRotor{iPatch} = {linspace(0, 1, q1), linspace(0, 1, q2)};
            end
            % Stator
            A_ges_stator = sum(obj.PatchAreaStator);
            for iPatch = 1:obj.NumberPatchesStator
                % lengths of sides L1: left+right, L2: top+bottom
                L1 = sum(msh_eval_boundary_side(obj.MultipatchMeshStator.msh_patch{iPatch}, 1).element_size) + ...
                   sum(msh_eval_boundary_side(obj.MultipatchMeshStator.msh_patch{iPatch}, 2).element_size); 
                L2 = sum(msh_eval_boundary_side(obj.MultipatchMeshStator.msh_patch{iPatch}, 3).element_size) + ...
                   sum(msh_eval_boundary_side(obj.MultipatchMeshStator.msh_patch{iPatch}, 4).element_size);
                num_quivers_i = ceil(obj.PatchAreaStator(iPatch)/A_ges_stator*obj.NumberQuiversStator);
                % Number quivers in x an y according to lengths of sides
                q1 = ceil((num_quivers_i* L2/L1)^0.5);
                q2 = ceil((num_quivers_i* L1/L2)^0.5);
                % Points with also half distance to boundary
                obj.ParamPlotPointsQuiverStator{iPatch} = {linspace(0, 1-1/q1, q1)+0.5/q1, linspace(0, 1-1/q2, q2)+0.5/q2};
                % Points also directly on boundary
                %obj.ParamPlotPointsQuiverStator{iPatch} = {linspace(0, 1, q1), linspace(0, 1, q2)};
            end
        end

        %% Plot the resulting B-Field
        function plotBResulting(obj)
            if obj.firstBplot == true
                obj.calcShapeValues();
                obj.firstBplot = false;
            end
            obj.calculateGeoPlotPoints();
            cols = 0.0:0.1:2;
            % Figure already exists, update current one
            if ishandle(obj.FigureBRes)
                % Rotor
                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    Ni_x_ui = sum(squeeze(obj.RotorShapeValues{iPatch}{2}(1, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    Ni_y_ui = sum(squeeze(obj.RotorShapeValues{iPatch}{2}(2, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    B_loc = (Ni_x_ui.^2 + Ni_y_ui.^2).^0.5;
                    set(obj.FigureBResBMagRotor{iPatch}, 'XData', obj.GeoPlotPointsRotor{iPatch}{1, 1}, 'YData', obj.GeoPlotPointsRotor{iPatch}{2, 1}, 'ZData', B_loc);
                end
                % Stator
                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    Ni_x_ui = sum(squeeze(obj.StatorShapeValues{iPatch}{2}(1, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    Ni_y_ui = sum(squeeze(obj.StatorShapeValues{iPatch}{2}(2, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    B_loc = (Ni_x_ui.^2 + Ni_y_ui.^2).^0.5;
                    set(obj.FigureBResBMagStator{iPatch}, 'XData', obj.GeoPlotPointsStator{iPatch}{1, 1}, 'YData', obj.GeoPlotPointsStator{iPatch}{2, 1}, 'ZData', B_loc);
                end
                % Boundaries and Interfaces
                % Rotor
                for iBnd = 1:numel(obj.BoundariesRotor)
                    patch = obj.BoundariesRotor(iBnd).patches;
                    side = obj.BoundariesRotor(iBnd).faces;
                    set(obj.FigureBResBoundariesRotor{iBnd}, 'XData', obj.GeoPlotPointsRotor{patch}{1, side+1}, 'YData', obj.GeoPlotPointsRotor{patch}{2, side+1});
                end
                for iInt = 1:numel(obj.InterfacesRotor)
                    patch = obj.InterfacesRotor(iInt).patch1;
                    side = obj.InterfacesRotor(iInt).side1;
                    set(obj.FigureBResInterfacesRotor{iInt}, 'XData', obj.GeoPlotPointsRotor{patch}{1, side+1}, 'YData', obj.GeoPlotPointsRotor{patch}{2, side+1});
                end
                % Stator
                for iBnd = 1:numel(obj.BoundariesStator)
                    patch = obj.BoundariesStator(iBnd).patches;
                    side = obj.BoundariesStator(iBnd).faces;
                    set(obj.FigureBResBoundariesStator{iBnd}, 'XData', obj.GeoPlotPointsStator{patch}{1, side+1}, 'YData', obj.GeoPlotPointsStator{patch}{2, side+1});
                end
                for iInt = 1:numel(obj.InterfacesStator)
                    patch = obj.InterfacesStator(iInt).patch1;
                    side = obj.InterfacesStator(iInt).side1;
                    set(obj.FigureBResInterfacesStator{iInt}, 'XData', obj.GeoPlotPointsStator{patch}{1, side+1}, 'YData', obj.GeoPlotPointsStator{patch}{2, side+1});
                end
                % Quivers
                Bx_vals = [];
                By_vals = [];
                X_vals = [];
                Y_vals =[];
                % Rotor quiver
                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    Bx = sum(obj.RotorShapeValuesQuiver{iPatch}{2}(2, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    By = -sum(obj.RotorShapeValuesQuiver{iPatch}{2}(1, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    x = obj.GeoPlotPointsQuiverRotor{iPatch}{1, 1};
                    y = obj.GeoPlotPointsQuiverRotor{iPatch}{2, 1};

                    X_vals = cat(2, X_vals, reshape(x, 1, []));
                    Y_vals = cat(2, Y_vals, reshape(y, 1, []));
                    Bx_vals = cat(2, Bx_vals, reshape(Bx, 1, []));
                    By_vals = cat(2, By_vals, reshape(By, 1, []));
                end
                % Rotate Rotor quivers
                Bxtemp = Bx_vals*cos(obj.RotationAngle) - By_vals*sin(obj.RotationAngle);
                Bytemp = By_vals*cos(obj.RotationAngle) + Bx_vals*sin(obj.RotationAngle);
                Bx_vals = Bxtemp;
                By_vals = Bytemp;
                % Stator quiver
                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    Bx = sum(obj.StatorShapeValuesQuiver{iPatch}{2}(2, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    By = -sum(obj.StatorShapeValuesQuiver{iPatch}{2}(1, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    x = obj.GeoPlotPointsQuiverStator{iPatch}{1, 1};
                    y = obj.GeoPlotPointsQuiverStator{iPatch}{2, 1};

                    X_vals = cat(2, X_vals, reshape(x, 1, []));
                    Y_vals = cat(2, Y_vals, reshape(y, 1, []));
                    Bx_vals = cat(2, Bx_vals, reshape(Bx, 1, []));
                    By_vals = cat(2, By_vals, reshape(By, 1, []));
                end
                set(obj.FigureBResQuiver, 'XData', X_vals, 'YData', Y_vals, 'UData', Bx_vals, 'VData', By_vals);
                set(obj.FigureBRes.CurrentAxes, "XLim", obj.XLimits)
                set(obj.FigureBRes.CurrentAxes, "YLim", obj.YLimits)
                drawnow();

            else %Create new Figure and plot 
                obj.FigureBRes = figure('Name', 'B Resulting Plot', 'NumberTitle', 'off', 'Position', [0 0 800 500]);
                hold on
                % Rotor
                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    Ni_x_ui = sum(squeeze(obj.RotorShapeValues{iPatch}{2}(1, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    Ni_y_ui = sum(squeeze(obj.RotorShapeValues{iPatch}{2}(2, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    B_loc = (Ni_x_ui.^2 + Ni_y_ui.^2).^0.5;
                    [~, obj.FigureBResBMagRotor{iPatch}] = contourf(obj.GeoPlotPointsRotor{iPatch}{1, 1}, obj.GeoPlotPointsRotor{iPatch}{2, 1}, B_loc, cols);
                    set(obj.FigureBResBMagRotor{iPatch}, 'LineColor', 'none');
                end
                % Stator
                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    Ni_x_ui = sum(squeeze(obj.StatorShapeValues{iPatch}{2}(1, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    Ni_y_ui = sum(squeeze(obj.StatorShapeValues{iPatch}{2}(2, :, :, :)).* reshape(u_loc, 1, 1, []), 3);
                    B_loc = (Ni_x_ui.^2 + Ni_y_ui.^2).^0.5;
                    [~, obj.FigureBResBMagStator{iPatch}] = contourf(obj.GeoPlotPointsStator{iPatch}{1, 1}, obj.GeoPlotPointsStator{iPatch}{2, 1}, B_loc, cols);
                    set(obj.FigureBResBMagStator{iPatch}, 'LineColor', 'none');
                end

                % Patch Boundaries
                for iBnd = 1:numel(obj.BoundariesRotor)
                    patch = obj.BoundariesRotor(iBnd).patches;
                    side = obj.BoundariesRotor(iBnd).faces;
                    obj.FigureBResBoundariesRotor{iBnd} = plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                for iInt = 1:numel(obj.InterfacesRotor)
                    patch = obj.InterfacesRotor(iInt).patch1;
                    side = obj.InterfacesRotor(iInt).side1;
                    obj.FigureBResInterfacesRotor{iInt} = plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                % Stator
                for iBnd = 1:numel(obj.BoundariesStator)
                    patch = obj.BoundariesStator(iBnd).patches;
                    side = obj.BoundariesStator(iBnd).faces;
                    obj.FigureBResBoundariesStator{iBnd} = plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                for iInt = 1:numel(obj.InterfacesStator)
                    patch = obj.InterfacesStator(iInt).patch1;
                    side = obj.InterfacesStator(iInt).side1;
                    obj.FigureBResInterfacesStator{iInt} = plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                % Quivers
                Bx_vals = [];
                By_vals = [];
                X_vals = [];
                Y_vals =[];
                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    Bx = sum(obj.RotorShapeValuesQuiver{iPatch}{2}(2, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    By = -sum(obj.RotorShapeValuesQuiver{iPatch}{2}(1, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    %
                    x = obj.GeoPlotPointsQuiverRotor{iPatch}{1, 1};
                    y = obj.GeoPlotPointsQuiverRotor{iPatch}{2, 1};

                    X_vals = cat(2, X_vals, reshape(x, 1, []));
                    Y_vals = cat(2, Y_vals, reshape(y, 1, []));
                    Bx_vals = cat(2, Bx_vals, reshape(Bx, 1, []));
                    By_vals = cat(2, By_vals, reshape(By, 1, []));
                end
                % Rotate Rotor quivers
                Bxtemp = Bx_vals*cos(obj.RotationAngle) - By_vals*sin(obj.RotationAngle);
                Bytemp = By_vals*cos(obj.RotationAngle) + Bx_vals*sin(obj.RotationAngle);
                Bx_vals = Bxtemp;
                By_vals = Bytemp;

                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    Bx = sum(obj.StatorShapeValuesQuiver{iPatch}{2}(2, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    By = -sum(obj.StatorShapeValuesQuiver{iPatch}{2}(1, :, :, :).* reshape(u_loc, 1, 1, 1, []), 4);
                    %
                    x = obj.GeoPlotPointsQuiverStator{iPatch}{1, 1};
                    y = obj.GeoPlotPointsQuiverStator{iPatch}{2, 1};

                    X_vals = cat(2, X_vals, reshape(x, 1, []));
                    Y_vals = cat(2, Y_vals, reshape(y, 1, []));
                    Bx_vals = cat(2, Bx_vals, reshape(Bx, 1, []));
                    By_vals = cat(2, By_vals, reshape(By, 1, []));
                end

                obj.FigureBResQuiver = quiver(X_vals, Y_vals, Bx_vals, By_vals, 'color', 'k');
                clims = [min(cols), max(cols)];
                caxis(clims);
                c = colorbar();
                c.TickLabelInterpreter = 'latex';
                c.Label.Interpreter = 'latex';
                c.Label.String = 'Magnetic flux density (T)';
                c.FontSize = 14;
                c.Label.FontSize = 16;
                view([0, 90]);
                axis equal;
                set(obj.FigureBRes.CurrentAxes, "XLim", obj.XLimits)
                set(obj.FigureBRes.CurrentAxes, "YLim", obj.YLimits)
            end
        end

        %% Plot the resulting Magnetic potential Lines
        function plotMagneticPotentialLines(obj)
            if obj.firstBplot == true
                obj.calcShapeValues();
                obj.firstBplot = false;
            end
            obj.calculateGeoPlotPoints();

            if ishandle(obj.FigurePotLines)
                % Update Geometry Surface and Boundaries
                % Rotor
                for iBnd = 1:numel(obj.BoundariesRotor)
                    patch = obj.BoundariesRotor(iBnd).patches;
                    side = obj.BoundariesRotor(iBnd).faces;
                    set(obj.FigurePotLinesBoundariesRotor{iBnd}, 'XData', obj.GeoPlotPointsRotor{patch}{1, side+1}, 'YData', obj.GeoPlotPointsRotor{patch}{2, side+1});
                end
                for iInt = 1:numel(obj.InterfacesRotor)
                    patch = obj.InterfacesRotor(iInt).patch1;
                    side = obj.InterfacesRotor(iInt).side1;
                    set(obj.FigurePotLinesInterfacesRotor{iInt}, 'XData', obj.GeoPlotPointsRotor{patch}{1, side+1}, 'YData', obj.GeoPlotPointsRotor{patch}{2, side+1});
                end
                % Stator
                for iBnd = 1:numel(obj.BoundariesStator)
                    patch = obj.BoundariesStator(iBnd).patches;
                    side = obj.BoundariesStator(iBnd).faces;
                    set(obj.FigurePotLinesBoundariesStator{iBnd}, 'XData', obj.GeoPlotPointsStator{patch}{1, side+1}, 'YData', obj.GeoPlotPointsStator{patch}{2, side+1});
                end
                for iInt = 1:numel(obj.InterfacesStator)
                    patch = obj.InterfacesStator(iInt).patch1;
                    side = obj.InterfacesStator(iInt).side1;
                    set(obj.FigurePotLinesInterfacesStator{iInt}, 'XData', obj.GeoPlotPointsStator{patch}{1, side+1}, 'YData', obj.GeoPlotPointsStator{patch}{2, side+1});
                end
                % Update Potential Lines
                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    A_loc = obj.RotorShapeValues{iPatch}{1}.* reshape(u_loc, 1, 1, []);
                    A_loc = abs(sum(A_loc, 3));
%                     B_loc = sum(obj.RotorShapeValues{iPatch}{2}.* reshape(u_loc, 1, 1, 1, []), 4);
%                     B_loc = squeeze((B_loc(1, :, :).^2 + B_loc(2, :, :).^2).^0.5);
                    set(obj.FigurePotLinesGeometryRotor{iPatch}, 'XData', obj.GeoPlotPointsRotor{iPatch}{1, 1}, 'YData', obj.GeoPlotPointsRotor{iPatch}{2, 1});
                    set(obj.FigurePotLinesPotentialLinesRotor{iPatch}, 'XData', obj.GeoPlotPointsRotor{iPatch}{1, 1}, 'YData', obj.GeoPlotPointsRotor{iPatch}{2, 1}, 'ZData', A_loc);
                    %surf(obj.GeoPlotPointsRotor{iPatch}{1, 1}, obj.GeoPlotPointsRotor{iPatch}{2, 1}, A_loc);
                end
                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    A_loc = obj.StatorShapeValues{iPatch}{1}.* reshape(u_loc, 1, 1, []);
                    A_loc = abs(sum(A_loc, 3));
%                     B_loc = sum(obj.StatorShapeValues{iPatch}{2}.* reshape(u_loc, 1, 1, 1, []), 4);
%                     B_loc = squeeze((B_loc(1, :, :).^2 + B_loc(2, :, :).^2).^0.5);
                    set(obj.FigurePotLinesGeometryStator{iPatch}, 'XData', obj.GeoPlotPointsStator{iPatch}{1, 1}, 'YData', obj.GeoPlotPointsStator{iPatch}{2, 1});
                    set(obj.FigurePotLinesPotentialLinesStator{iPatch}, 'XData', obj.GeoPlotPointsStator{iPatch}{1, 1}, 'YData', obj.GeoPlotPointsStator{iPatch}{2, 1}, 'ZData', A_loc);
                    %surf(obj.GeoPlotPointsStator{iPatch}{1, 1}, obj.GeoPlotPointsStator{iPatch}{2, 1}, A_loc);
                end
                drawnow();
            else
                % Create new plot for Geometry
                obj.FigurePotLines = figure('Name', 'Potential Line Plot', 'NumberTitle', 'off', 'Position', [0 0 800 500]);
                hold on
                axis equal
                % Draw Geometry Rotor
                for iMat = 1:numel(obj.Materials)
                    for patch = obj.Materials(iMat).PatchesRotor
                        obj.FigurePotLinesGeometryRotor{patch}=surf(obj.GeoPlotPointsRotor{patch}{1, 1}, obj.GeoPlotPointsRotor{patch}{2, 1}, zeros(size(obj.GeoPlotPointsRotor{patch}{1, 1})), ...
                            "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.PlotFaceAlpha);
                    end
                    for patch = obj.Materials(iMat).PatchesStator
                        obj.FigurePotLinesGeometryStator{patch}=surf(obj.GeoPlotPointsStator{patch}{1, 1}, obj.GeoPlotPointsStator{patch}{2, 1}, zeros(size(obj.GeoPlotPointsStator{patch}{1, 1})), ...
                            "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.PlotFaceAlpha);
                    end
                end
                % Boundaries
                for iBnd = 1:numel(obj.BoundariesRotor)
                    patch = obj.BoundariesRotor(iBnd).patches;
                    side = obj.BoundariesRotor(iBnd).faces;
                    obj.FigurePotLinesBoundariesRotor{iBnd} = plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                for iInt = 1:numel(obj.InterfacesRotor)
                    patch = obj.InterfacesRotor(iInt).patch1;
                    side = obj.InterfacesRotor(iInt).side1;
                    obj.FigurePotLinesInterfacesRotor{iInt} = plot(obj.GeoPlotPointsRotor{patch}{1, side+1}, obj.GeoPlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                % Stator
                for iBnd = 1:numel(obj.BoundariesStator)
                    patch = obj.BoundariesStator(iBnd).patches;
                    side = obj.BoundariesStator(iBnd).faces;
                    obj.FigurePotLinesBoundariesStator{iBnd} = plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end
                for iInt = 1:numel(obj.InterfacesStator)
                    patch = obj.InterfacesStator(iInt).patch1;
                    side = obj.InterfacesStator(iInt).side1;
                    obj.FigurePotLinesInterfacesStator{iInt} = plot(obj.GeoPlotPointsStator{patch}{1, side+1}, obj.GeoPlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                end

                % Calculate Potential Lines
                min_val = Inf;
                max_val = -Inf;
                % get min and max values for consistency over all patches
                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    A_loc = obj.RotorShapeValues{iPatch}{1}.* reshape(u_loc, 1, 1, []);
                    A_loc = abs(sum(A_loc, 3));
%                     B_loc = sum(obj.RotorShapeValues{iPatch}{2}.* reshape(u_loc, 1, 1, 1, []), 4);
%                     B_loc = squeeze((B_loc(1, :, :).^2 + B_loc(2, :, :).^2).^0.5);
                    min_val = min(min(A_loc, [], [1, 2]), min_val);
                    max_val = max(max(A_loc, [], [1, 2]), max_val);
                end
                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    A_loc = obj.StatorShapeValues{iPatch}{1}.* reshape(u_loc, 1, 1, []);
                    A_loc = abs(sum(A_loc, 3));
%                     B_loc = sum(obj.StatorShapeValues{iPatch}{2}.* reshape(u_loc, 1, 1, 1, []), 4);
%                     B_loc = squeeze((B_loc(1, :, :).^2 + B_loc(2, :, :).^2).^0.5);
                    min_val = min(min(A_loc, [], [1, 2]), min_val);
                    max_val = max(max(A_loc, [], [1, 2]), max_val);
                end
                obj.FigurePotLinesMin = min_val;
                obj.FigurePotLinesMax = max_val;   

                for iPatch = 1:obj.NumberPatchesRotor
                    u_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                    A_loc = obj.RotorShapeValues{iPatch}{1}.* reshape(u_loc, 1, 1, []);
                    A_loc = abs(sum(A_loc, 3));
%                     B_loc = sum(obj.RotorShapeValues{iPatch}{2}.* reshape(u_loc, 1, 1, 1, []), 4);
%                     B_loc = squeeze((B_loc(1, :, :).^2 + B_loc(2, :, :).^2).^0.5);
                    [~, obj.FigurePotLinesPotentialLinesRotor{iPatch}] = contour(obj.GeoPlotPointsRotor{iPatch}{1, 1}, obj.GeoPlotPointsRotor{iPatch}{2, 1}, A_loc, linspace(min_val, max_val, obj.NumPotLines));%, LineColor="black");
                end
                for iPatch = 1:obj.NumberPatchesStator
                    u_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                    A_loc = obj.StatorShapeValues{iPatch}{1}.* reshape(u_loc, 1, 1, []);
                    A_loc = abs(sum(A_loc, 3));
%                     B_loc = sum(obj.StatorShapeValues{iPatch}{2}.* reshape(u_loc, 1, 1, 1, []), 4);
%                     B_loc = squeeze((B_loc(1, :, :).^2 + B_loc(2, :, :).^2).^0.5);
                    [~, obj.FigurePotLinesPotentialLinesStator{iPatch}] = contour(obj.GeoPlotPointsStator{iPatch}{1, 1}, obj.GeoPlotPointsStator{iPatch}{2, 1}, A_loc, linspace(min_val, max_val, obj.NumPotLines));%, LineColor="black");
                 end
            
                %shading interp
                colormap jet
                view([0, 90]);
                axis equal;
                set(obj.FigurePotLines.CurrentAxes, "XLim", obj.XLimits)
                set(obj.FigurePotLines.CurrentAxes, "YLim", obj.YLimits)
            end
        end

        %% Plot Gap Evaluation Points
        function plotEvalPoints(obj, draw_text)
            if nargin == 1
               draw_text = false; 
            end
            for iPatch = 1:numel(obj.EvalPatchesRotor)
                x_physical = obj.MultipatchMeshRotor.msh_patch{obj.EvalPatchesRotor(iPatch)}.map(obj.EvalPointsRotor(:, iPatch));
                scatter(x_physical(1), x_physical(2), "black", "filled");
                if draw_text
                    text(x_physical(1)+0.001, x_physical(2), num2str(iPatch));
                end
            end
            for iPatch = 1:numel(obj.EvalPatchesStator)
                x_physical = obj.MultipatchMeshStator.msh_patch{obj.EvalPatchesStator(iPatch)}.map(obj.EvalPointsStator(:, iPatch));
                scatter(x_physical(1), x_physical(2), "black", "filled");
                if draw_text
                    text(x_physical(1)+0.001, x_physical(2), num2str(iPatch));
                end
            end
        end

        function plotRemanenceQuiver(obj)
            xvals = [];
            yvals = [];
            Bx = [];
            By = [];
            for remPatch = find(~cellfun(@isempty, obj.PlotRemanence))
                xvals = [xvals, reshape(obj.GeoPlotPointsQuiverRotor{remPatch}{1, 1}, 1, [])];
                yvals = [yvals, reshape(obj.GeoPlotPointsQuiverRotor{remPatch}{2, 1}, 1, [])];
                Brx = reshape(obj.PlotRemanence{remPatch}(:, 1), 1, []);
                Bry = reshape(obj.PlotRemanence{remPatch}(:, 2), 1, []);
                Bx = [Bx, Brx*cos(obj.RotationAngle) - Bry*sin(obj.RotationAngle)];
                By = [By, Brx*sin(obj.RotationAngle) + Bry*cos(obj.RotationAngle)];
            end
            quiver(xvals, yvals, Bx, By, 3e-1, 'color', 'k', 'LineWidth', 1.5);
        end
        
        %% Export Magnetic Potential to ParaView
        function exportPotentialParaView(obj, filename, angleIndex, angleList)
            if nargin <= 1
                filename = 'paraview/MagMotor';
            end
            if nargin <= 2
                angleIndex = 1;
                angleList = 1;
            end

            str1 = cat (2, '<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
            str2 = cat (2, '<DataSet timestep="%d" group="" part="%d" file="%s.vts"/> \n');
            str3 = cat (2, '</Collection>\n', '</VTKFile> \n');

            if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
                pvd_filename = cat (2, filename, '.pvd');
            else
                pvd_filename = filename;
                filename = filename (1:end-4);
            end

            if angleIndex == 1
                fid = fopen (pvd_filename, 'w');
                fprintf (fid, str1);
            else
                fid = fopen (pvd_filename, 'a');
            end
            
            if (fid < 0)
                error ('exportSolutionParaView: could not open file %s', pvd_filename);
            end

            ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
            if (isempty (ind)); ind = 1; end
            obj.calculateGeoPlotPoints;
            % Rotor
            for iPatch = 1:obj.NumberPatchesRotor
                filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(angleIndex-1), '_', num2str (iPatch));
                filename_patch = cat (2, filename, '_', num2str(angleIndex-1), '_', num2str (iPatch));
                vts_pts = [];
                %%
                vts_pts(1, :, :) = obj.GeoPlotPointsRotor{iPatch}{1, 1};
                vts_pts(2, :, :) = obj.GeoPlotPointsRotor{iPatch}{2, 1};
                M_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                M_pts = sum(obj.RotorShapeValues{iPatch}{1}(:, :, :).* reshape(M_loc, 1, 1, []), 3);

                msh_to_vtk (vts_pts, M_pts, filename_patch, 'A_z');

                fprintf (fid, str2, angleList(angleIndex), iPatch, filename_patch_without_path);

                
            end
            % Stator
            for iPatch = 1:obj.NumberPatchesStator
                filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(angleIndex-1), '_', num2str (iPatch+obj.NumberPatchesRotor));
                filename_patch = cat (2, filename, '_', num2str(angleIndex-1), '_', num2str (iPatch+obj.NumberPatchesRotor));
                vts_pts = [];
                vts_pts(1, :, :) = obj.GeoPlotPointsStator{iPatch}{1, 1};
                vts_pts(2, :, :) = obj.GeoPlotPointsStator{iPatch}{2, 1};
                M_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                M_pts = sum(obj.StatorShapeValues{iPatch}{1}(:, :, :).* reshape(M_loc, 1, 1, []), 3);
                msh_to_vtk (vts_pts, M_pts, filename_patch, 'A_z');
                fprintf (fid, str2, angleList(angleIndex), iPatch+obj.NumberPatchesRotor, filename_patch_without_path);
            end
            if angleIndex == numel(angleList)
                fprintf (fid, str3);
            end
            fclose (fid);

        end

        %% Export Magnetic B Field to ParaView
        function exportBFieldParaView(obj, filename, currentIndex, dataList)
            if nargin <= 1
                filename = 'paraview/BFieldMotor';
            end
            if nargin <= 2
                currentIndex = 1;
                dataList = 1;
            end

            str1 = cat (2, '<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
            str2 = cat (2, '<DataSet timestep="%d" group="" part="%d" file="%s.vts"/> \n');
            str3 = cat (2, '</Collection>\n', '</VTKFile> \n');

            if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
                pvd_filename = cat (2, filename, '.pvd');
            else
                pvd_filename = filename;
                filename = filename (1:end-4);
            end

            if currentIndex == 1
                fid = fopen (pvd_filename, 'w');
                fprintf (fid, str1);
            else
                fid = fopen (pvd_filename, 'a');
            end
            
            if (fid < 0)
                error ('exportSolutionParaView: could not open file %s', pvd_filename);
            end

            ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
            if (isempty (ind)); ind = 1; end
            obj.calculateGeoPlotPoints;
            % Rotor
            for iPatch = 1:obj.NumberPatchesRotor
                filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(currentIndex-1), '_', num2str (iPatch));
                filename_patch = cat (2, filename, '_', num2str(currentIndex-1), '_', num2str (iPatch));
                vts_pts = [];
                %%
                vts_pts(1, :, :) = obj.GeoPlotPointsQuiverRotor{iPatch}{1, 1};
                vts_pts(2, :, :) = obj.GeoPlotPointsQuiverRotor{iPatch}{2, 1};
                M_loc = obj.MagneticPotentialRotor(obj.MultipatchSpaceRotor.gnum{iPatch});
                Bx = sum(obj.RotorShapeValuesQuiver{iPatch}{2}(2, :, :, :).* reshape(M_loc, 1, 1, 1, []), 4);
                By = -sum(obj.RotorShapeValuesQuiver{iPatch}{2}(1, :, :, :).* reshape(M_loc, 1, 1, 1, []), 4);
                
                Bxtemp = Bx*cos(obj.RotationAngle) - By*sin(obj.RotationAngle);
                Bytemp = By*cos(obj.RotationAngle) + Bx*sin(obj.RotationAngle);
                B = [Bxtemp;Bytemp];

                msh_to_vtk_vec (vts_pts, B, filename_patch, 'B');
                fprintf (fid, str2, dataList(currentIndex), iPatch, filename_patch_without_path);

                
            end
            
            % Stator
            for iPatch = 1:obj.NumberPatchesStator
                filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(currentIndex-1), '_', num2str (iPatch+obj.NumberPatchesRotor));
                filename_patch = cat (2, filename, '_', num2str(currentIndex-1), '_', num2str (iPatch+obj.NumberPatchesRotor));
                vts_pts = [];
                vts_pts(1, :, :) = obj.GeoPlotPointsQuiverStator{iPatch}{1, 1};
                vts_pts(2, :, :) = obj.GeoPlotPointsQuiverStator{iPatch}{2, 1};
                M_loc = obj.MagneticPotentialStator(obj.MultipatchSpaceStator.gnum{iPatch});
                Bx = sum(obj.StatorShapeValuesQuiver{iPatch}{2}(2, :, :, :).* reshape(M_loc, 1, 1, 1, []), 4);
                By = -sum(obj.StatorShapeValuesQuiver{iPatch}{2}(1, :, :, :).* reshape(M_loc, 1, 1, 1, []), 4);
                B = [Bx;By];

                msh_to_vtk_vec (vts_pts, B, filename_patch, 'B');
                fprintf (fid, str2, dataList(currentIndex), iPatch+obj.NumberPatchesRotor, filename_patch_without_path);
            end
            
            if currentIndex == numel(dataList)
                fprintf (fid, str3);
            end
            fclose (fid);

        end

        function plotCouplingValues(obj)
            figure(12)
            clf
            sinValues = obj.CouplingValues((mod(obj.CouplingIndices, 2) == 1));
            cosValues = obj.CouplingValues((mod(obj.CouplingIndices, 2) == 0));
            xlims = [min(obj.HarmonicsAll)-1, max(obj.HarmonicsAll)+1];
            subplot(2,1,1)
            bar(obj.HarmonicsSin, sinValues);
            xticks(obj.HarmonicsSin);
            xlim(xlims);
            title("Sin Values");
            subplot(2,1,2)
            bar(obj.HarmonicsCos, cosValues);
            xticks(obj.HarmonicsCos);
            xlim(xlims);
            title("Cos Values");
        end

        %% Param plot points
        function calculateParamPlotPoints(obj)
            disp("Calculating Parametric Plot Points...")
            for iPatch = 1:obj.NumberPatchesRotor
                pointsX = unique(kntrefine(obj.GeometryRotor(iPatch).nurbs.knots{1}, obj.PlotResolutionX, obj.GeometryRotor(iPatch).nurbs.order(1), obj.GeometryRotor(iPatch).nurbs.order(1)-1));
                pointsY = unique(kntrefine(obj.GeometryRotor(iPatch).nurbs.knots{2}, obj.PlotResolutionY, obj.GeometryRotor(iPatch).nurbs.order(2), obj.GeometryRotor(iPatch).nurbs.order(2)-1));
                obj.ParamPlotPointsRotor{iPatch} = {pointsX, pointsY};
            end
            for iPatch = 1:obj.NumberPatchesStator
                pointsX = unique(kntrefine(obj.GeometryStator(iPatch).nurbs.knots{1}, obj.PlotResolutionX, obj.GeometryStator(iPatch).nurbs.order(1), obj.GeometryStator(iPatch).nurbs.order(1)-1));
                pointsY = unique(kntrefine(obj.GeometryStator(iPatch).nurbs.knots{2}, obj.PlotResolutionY, obj.GeometryStator(iPatch).nurbs.order(2), obj.GeometryStator(iPatch).nurbs.order(2)-1));
                obj.ParamPlotPointsStator{iPatch} = {pointsX, pointsY};
            end
        end

        %% Calculate geometry points
        function calculateGeoPlotPoints(obj)
%             disp("Calculating Physical Plot Points...")
            for iPatch = 1:obj.NumberPatchesRotor
                obj.GeoPlotPointsRotor{iPatch} = obj.getGeoPlotPoints(obj.GeometryRotor(iPatch), obj.ParamPlotPointsRotor{iPatch}, obj.RotationAngle);
                obj.GeoPlotPointsQuiverRotor{iPatch} = obj.getGeoPlotPoints(obj.GeometryRotor(iPatch), obj.ParamPlotPointsQuiverRotor{iPatch}, obj.RotationAngle);
            end
            for iPatch = 1:obj.NumberPatchesStator
                obj.GeoPlotPointsStator{iPatch} = obj.getGeoPlotPoints(obj.GeometryStator(iPatch), obj.ParamPlotPointsStator{iPatch}, 0);
                obj.GeoPlotPointsQuiverStator{iPatch} = obj.getGeoPlotPoints(obj.GeometryStator(iPatch), obj.ParamPlotPointsQuiverStator{iPatch}, 0);
            end
        end

        function pnts = getGeoPlotPoints(obj, geo, pts, rot)
            xpts = pts{1};
            ypts = pts{2};
        
            F_geo = nrbeval(geo.nurbs, {xpts, ypts});
            X = squeeze(F_geo(1, :, :))*cos(rot) - squeeze(F_geo(2, :, :))*sin(rot);
            Y = squeeze(F_geo(2, :, :))*cos(rot) + squeeze(F_geo(1, :, :))*sin(rot);
            pnts{1, 1} = X; 
            pnts{2, 1} = Y; 
        
            bnd = nrbeval(geo.boundary(1).nurbs, ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 2} = X;
            pnts{2, 2} = Y;
        
            bnd = nrbeval(geo.boundary(2).nurbs, ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 3} = X;
            pnts{2, 3} = Y;
        
            bnd = nrbeval(geo.boundary(3).nurbs, xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 4} = X;
            pnts{2, 4} = Y;
        
            bnd = nrbeval(geo.boundary(4).nurbs, xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 5} = X;
            pnts{2, 5} = Y;  
        end

    end
end
