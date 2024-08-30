classdef MotorOptimizationFull < MotorSimulation

    properties
        % Optimization
        NumOptParameters=0;
        NumOptControlPointsRotor=0;
        NumOptControlPointsStator=0;
        StartOptTime, SaveDirectory, OptimizationFigure;
        srfplotRotor, lplotRotor, intplotRotor;
        srfplotStator, lplotStator, intplotStator;
        Options;
        RotorGeometryFunction, StatorGeometryFunction;
        ConstraintFunction;
        OptimizationParameters;
        OptimizationControlPointsRotor, OptimizationControlPointsStator;
        UpperBounds, LowerBounds;
        LastX, LastT, LastdT;
        OptimizationAngles;
        OptimizationType="linear";
        MinOptiTorque;
        OptimizationHistory=struct('iteration', {}, 'funccount', {}, 'fval', {}, 'fvalUsed', {}, 'Tmean', {}, 'Tsigma', {}, 'Cost', {}, 'Currents',  {}, 'Angles', {}, 'Torques', {}, 'Parameters', {}, 'TimeElapsed', {}, 'Power', {}, 'Phases', {}, 'Mult', {});
        MultCost, MultTsig;
        MultPower;
        OptimizationMethod = "interior-point";
        TimeStart;
        Jinit = 3.2e+06;  % Current density in phase 1 of initial design.
%         OptimizationPhases
    end

    methods
        %% constructor, set material properties
        function obj = MotorOptimizationFull()
            obj@MotorSimulation();
            obj.resetOptimizationParameters();
            obj.resetOptimizationControlPoints();
            obj.TimeStart = tic;
        end

        function resetOptimizationParameters(obj)
            obj.NumOptParameters = 0;
            obj.OptimizationParameters=struct('Name', {}, 'InitVal', {}, 'MinVal', {}, 'MaxVal', {}, 'OptVal', {}, 'MagnetPatches', {}, 'MagnetRelations', {});
        end

        function addOptimizationParameter(obj, name, initVal, minVal, maxVal, MagnetPatches, MagnetRelations)
            obj.NumOptParameters = obj.NumOptParameters +1;
            obj.OptimizationParameters(obj.NumOptParameters).Name = name;
            obj.OptimizationParameters(obj.NumOptParameters).InitVal = initVal;
            obj.OptimizationParameters(obj.NumOptParameters).MinVal = minVal;
            obj.OptimizationParameters(obj.NumOptParameters).MaxVal = maxVal;
            if nargin > 5
                % Patches that contain a Permanent magnet
                obj.OptimizationParameters(obj.NumOptParameters).MagnetPatches = MagnetPatches;
                % Function for each patch on how to calculate the magnet
                % angle from the optimization parameter
                obj.OptimizationParameters(obj.NumOptParameters).MagnetRelations = MagnetRelations;
            end
            % Update all Min and Max Values
            obj.LowerBounds = [vertcat(obj.OptimizationParameters.MinVal); vertcat(obj.OptimizationControlPointsRotor.MinOffset); vertcat(obj.OptimizationControlPointsStator.MinOffset)];
            obj.UpperBounds = [vertcat(obj.OptimizationParameters.MaxVal); vertcat(obj.OptimizationControlPointsRotor.MaxOffset); vertcat(obj.OptimizationControlPointsStator.MaxOffset)]; 
        end

        function resetOptimizationControlPoints(obj)
            obj.NumOptControlPointsRotor = 0;
            obj.NumOptControlPointsStator = 0;
            obj.OptimizationControlPointsRotor=struct('Number', {}, 'Angle', {}, 'InitOffset', {}, 'MinOffset', {}, 'MaxOffset', {});
            obj.OptimizationControlPointsStator=struct('Number', {}, 'Angle', {}, 'InitOffset', {}, 'MinOffset', {}, 'MaxOffset', {});
        end
  
        function addOptimizationControlPointRotor(obj, pointNr, initOffset, minOffset, maxOffset)
            obj.NumOptControlPointsRotor = obj.NumOptControlPointsRotor +1;
            obj.OptimizationControlPointsRotor(obj.NumOptControlPointsRotor).Number = pointNr;
            obj.OptimizationControlPointsRotor(obj.NumOptControlPointsRotor).Angle = atan2(obj.ControlPointsRotor(pointNr, 2), obj.ControlPointsRotor(pointNr, 1));
            obj.OptimizationControlPointsRotor(obj.NumOptControlPointsRotor).InitOffset = initOffset;
            % Min and Max offsets of points
            obj.OptimizationControlPointsRotor(obj.NumOptControlPointsRotor).MinOffset = minOffset;
            obj.OptimizationControlPointsRotor(obj.NumOptControlPointsRotor).MaxOffset = maxOffset;
            % Update all Min and Max Values
            obj.LowerBounds = [vertcat(obj.OptimizationParameters.MinVal); vertcat(obj.OptimizationControlPointsRotor.MinOffset); vertcat(obj.OptimizationControlPointsStator.MinOffset)];
            obj.UpperBounds = [vertcat(obj.OptimizationParameters.MaxVal); vertcat(obj.OptimizationControlPointsRotor.MaxOffset); vertcat(obj.OptimizationControlPointsStator.MaxOffset)]; 
        end

        function addOptimizationControlPointStator(obj, pointNr, initOffset, minOffset, maxOffset)
            obj.NumOptControlPointsStator = obj.NumOptControlPointsStator +1;
            obj.OptimizationControlPointsStator(obj.NumOptControlPointsStator).Number = pointNr;
            obj.OptimizationControlPointsStator(obj.NumOptControlPointsStator).Angle = atan2(obj.ControlPointsStator(pointNr, 2), obj.ControlPointsStator(pointNr, 1));
            obj.OptimizationControlPointsStator(obj.NumOptControlPointsStator).InitOffset = initOffset;
            % Min and Max offsets of points
            obj.OptimizationControlPointsStator(obj.NumOptControlPointsStator).MinOffset = minOffset;
            obj.OptimizationControlPointsStator(obj.NumOptControlPointsStator).MaxOffset = maxOffset;
            % Update all Min and Max Values
            obj.LowerBounds = [vertcat(obj.OptimizationParameters.MinVal); vertcat(obj.OptimizationControlPointsRotor.MinOffset); vertcat(obj.OptimizationControlPointsStator.MinOffset)];
            obj.UpperBounds = [vertcat(obj.OptimizationParameters.MaxVal); vertcat(obj.OptimizationControlPointsRotor.MaxOffset); vertcat(obj.OptimizationControlPointsStator.MaxOffset)]; 
        end
        % Update Phases for new geometry
        function optimizationGeneratePhases(obj, spaceStator, spaceStatorEval, spaceStatorGeo, spaceStatorGeoEval, meshStator, meshStatorEval, scaleR)
            for iPhase = 1:numel(obj.Phases)
                obj.Phases(iPhase).Area = op_Omega_mp_eval(meshStatorEval, obj.Phases(iPhase).Patches);
                obj.Phases(iPhase).RHS = op_f_v_mp_eval(spaceStator, spaceStatorEval, meshStatorEval, obj.Phases(iPhase).Patches);
                obj.Phases(iPhase).J = obj.Jinit * obj.Phases(1).Area / obj.Phases(iPhase).Area * scaleR;
                obj.Phases(iPhase).Current = obj.Phases(iPhase).J*obj.Phases(iPhase).Area/obj.NumWindings;
                obj.Phases(iPhase).dbdC = op_D_DC_f_v_mp_eval(spaceStator, spaceStatorEval, spaceStatorGeo, spaceStatorGeoEval, meshStator, meshStatorEval, 1, obj.Phases(iPhase).Patches);
                obj.Phases(iPhase).dAdC = op_D_DC_Omega_mp_eval(spaceStatorGeo, spaceStatorGeoEval, meshStator, meshStatorEval, obj.Phases(iPhase).Patches);
            end
            obj.ApplicationCurrent = obj.Phases(1).Current; % Current should be constant for all phases!
        end
   
        function [fopt, gradf] = optimizationFunctionParameterShape(obj, x)
            % Rescale to real values
            x = x.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;
            % Initial parameters
            opts = obj.Options;
            opts.draw_geometry = false;
            for iParam  = 1:numel (obj.OptimizationParameters)
              opts.(obj.OptimizationParameters(iParam).Name) = x(iParam);
            end
            % Surface and control points of initial geometry
            % Rotor
            [srfInitRotor, ~, ~] = obj.RotorGeometryFunction(opts);
            ControlPointsInitRotor = sparse(size(obj.ControlPointsRotor, 1), size(obj.ControlPointsRotor, 2));
            for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                ControlPointsInitRotor(ind_loc, 1) = reshape(srfInitRotor(i).coefs(1, :, :, :)./srfInitRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitRotor(ind_loc, 2) = reshape(srfInitRotor(i).coefs(2, :, :, :)./srfInitRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitRotor(ind_loc, 3) = reshape(srfInitRotor(i).coefs(3, :, :, :)./srfInitRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitRotor(ind_loc, 4) = reshape(srfInitRotor(i).coefs(4, :, :, :), [], 1);
            end
            % Stator
            [srfInitStator, ~, ~] = obj.StatorGeometryFunction(opts);
            ControlPointsInitStator = sparse(size(obj.ControlPointsStator, 1), size(obj.ControlPointsStator, 2));
            for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                ControlPointsInitStator(ind_loc, 1) = reshape(srfInitStator(i).coefs(1, :, :, :)./srfInitStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitStator(ind_loc, 2) = reshape(srfInitStator(i).coefs(2, :, :, :)./srfInitStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitStator(ind_loc, 3) = reshape(srfInitStator(i).coefs(3, :, :, :)./srfInitStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitStator(ind_loc, 4) = reshape(srfInitStator(i).coefs(4, :, :, :), [], 1);
            end
            
            % Apply control point offsets
            ControlPointOffsetsRotor = x(obj.NumOptParameters+(1:obj.NumOptControlPointsRotor));
            for iCtrPt = 1:numel(ControlPointOffsetsRotor)
                valsx = ControlPointOffsetsRotor(iCtrPt).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                valsy = ControlPointOffsetsRotor(iCtrPt).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                ControlPointsInitRotor(numbers, 1) = ControlPointsInitRotor(numbers, 1) + valsx;
                ControlPointsInitRotor(numbers, 2) = ControlPointsInitRotor(numbers, 2) + valsy;
            end
            % Stator one
            ControlPointOffsetsStator = x(obj.NumOptParameters+obj.NumOptControlPointsRotor+(1:obj.NumOptControlPointsStator));
            for iCtrPt = 1:numel(ControlPointOffsetsStator)
                valsx = ControlPointOffsetsStator(iCtrPt).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                valsy = ControlPointOffsetsStator(iCtrPt).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                ControlPointsInitStator(numbers, 1) = ControlPointsInitStator(numbers, 1) + valsx;
                ControlPointsInitStator(numbers, 2) = ControlPointsInitStator(numbers, 2) + valsy;
            end
            % Update initial surface
            % Rotor
            for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                srfInitRotor(i).coefs(1, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 1).*ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
                srfInitRotor(i).coefs(2, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 2).*ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
                srfInitRotor(i).coefs(3, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 3).*ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
                srfInitRotor(i).coefs(4, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
            end
            % Stator
            for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                srfInitStator(i).coefs(1, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 1).*ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
                srfInitStator(i).coefs(2, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 2).*ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
                srfInitStator(i).coefs(3, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 3).*ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
                srfInitStator(i).coefs(4, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
            end
            
            % Load modified geometries: Rotor
            [GeometryRotor1, ~, ~, ~, ~] = mp_geo_load(srfInitRotor);
            meshes = cell (1, obj.NumberPatchesRotor);
            spaces  = cell (1, obj.NumberPatchesRotor);
            spaces_geo = cell (1, obj.NumberPatchesRotor);
            % create space partition and form functions
            for iptc = 1:obj.NumberPatchesRotor
              [knots, zeta] = kntrefine (GeometryRotor1(iptc).nurbs.knots, obj.SubdivisionsPatchesRotor-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
              rule      = msh_gauss_nodes (obj.DegreeQuadrature);
              [qn, qw]  = msh_set_quad_nodes (zeta, rule);
              meshes{iptc} = msh_cartesian (zeta, qn, qw, GeometryRotor1(iptc));
              spaces{iptc}  = sp_bspline (knots, obj.FormFunctionDegree, meshes{iptc});
              spaces_geo{iptc} = sp_nurbs(GeometryRotor1(iptc).nurbs, meshes{iptc});
            end
            MultipatchMeshRotor1 = msh_multipatch(meshes, obj.BoundariesRotor);
            MultipatchSpaceRotor1  = sp_multipatch(spaces, MultipatchMeshRotor1, obj.InterfacesRotor, obj.BoundaryInterfacesRotor);
            MultipatchSpaceRotorGeo1 = sp_multipatch(spaces_geo, MultipatchMeshRotor1, obj.InterfacesRotor, obj.BoundaryInterfacesRotor);
            % Precalculate Spaces and Meshes
            mshRotorEval = cell (1, obj.NumberPatchesRotor);
            spaceRotorEval  = cell (1, obj.NumberPatchesRotor);
            spaceRotorEvalGeo = cell (1, obj.NumberPatchesRotor);
            for iPatch = 1:obj.NumberPatchesRotor
                mshRotorEval{iPatch} = msh_precompute(MultipatchMeshRotor1.msh_patch{iPatch});
                spaceRotorEval{iPatch}  = sp_precompute(MultipatchSpaceRotor1.sp_patch{iPatch}, mshRotorEval{iPatch}, 'gradient', true);
                spaceRotorEvalGeo{iPatch}  = sp_precompute_param(MultipatchSpaceRotorGeo1.sp_patch{iPatch}, mshRotorEval{iPatch}, 'value', true, 'gradient', true);
            end

            [GeometryStator1, ~, ~, ~, ~] = mp_geo_load(srfInitStator);
            meshes = cell (1, obj.NumberPatchesStator);
            spaces  = cell (1, obj.NumberPatchesStator);
            spaces_geo = cell (1, obj.NumberPatchesStator);
            % create space partition and form functions
            for iptc = 1:obj.NumberPatchesStator
              [knots, zeta] = kntrefine (GeometryStator1(iptc).nurbs.knots, obj.SubdivisionsPatchesStator-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
              rule      = msh_gauss_nodes (obj.DegreeQuadrature);
              [qn, qw]  = msh_set_quad_nodes (zeta, rule);
              meshes{iptc} = msh_cartesian (zeta, qn, qw, GeometryStator1(iptc));
              spaces{iptc}  = sp_bspline (knots, obj.FormFunctionDegree, meshes{iptc});
              spaces_geo{iptc} = sp_nurbs(GeometryStator1(iptc).nurbs, meshes{iptc});
            end
            MultipatchMeshStator1 = msh_multipatch(meshes, obj.BoundariesStator);
            MultipatchSpaceStator1  = sp_multipatch(spaces, MultipatchMeshStator1, obj.InterfacesStator, obj.BoundaryInterfacesStator);
            MultipatchSpaceStatorGeo1 = sp_multipatch(spaces_geo, MultipatchMeshStator1, obj.InterfacesStator, obj.BoundaryInterfacesStator);
            % Precalculate Spaces and Meshes
            mshStatorEval = cell (1, obj.NumberPatchesStator);
            spaceStatorEval  = cell (1, obj.NumberPatchesStator);
            spaceStatorEvalGeo = cell (1, obj.NumberPatchesStator);
            for iPatch = 1:obj.NumberPatchesStator
                mshStatorEval{iPatch} = msh_precompute(MultipatchMeshStator1.msh_patch{iPatch});
                spaceStatorEval{iPatch}  = sp_precompute(MultipatchSpaceStator1.sp_patch{iPatch}, mshStatorEval{iPatch}, 'gradient', true);
                spaceStatorEvalGeo{iPatch}  = sp_precompute_param(MultipatchSpaceStatorGeo1.sp_patch{iPatch}, mshStatorEval{iPatch}, 'value', true, 'gradient', true);
            end

            ScaleR = 1;
            if isfield(opts, 'ScaleR')
                ScaleR = opts.ScaleR;
            end

            % Build RHS for rotor with updated magnetization depending on parameters that change magnetization direction
            RHS1 = sparse(MultipatchSpaceRotor1.ndof, 1);
            for iRem = 1:numel(obj.Remanence)
                remParams = zeros(size(obj.Remanence(iRem).Description));
                for iParam  = 1:numel(obj.Remanence(iRem).Description)
                    remParams(iParam) = opts.(obj.Remanence(iRem).Description(iParam));
                end
                alpha = obj.Remanence(iRem).Relation(remParams);
                if isempty(alpha)
                    alpha = obj.Remanence(iRem).InitAngle;
                end
                patch = obj.Remanence(iRem).Patch;
                Hc = obj.Remanence(iRem).Material.Br / obj.Remanence(iRem).Material.getMuLinear();
                RHS1 = RHS1 + Hc*op_gradv_n_bot_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, mshRotorEval, alpha, patch);
            end
            obj.RHSRotor = RHS1;

            % Reconstruct RHS of stator for updated geometry
            obj.optimizationGeneratePhases(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, ScaleR);

            if nargout == 2
                % Numerical derivatives dC/dP of parameter iOpt
                dCrtdP = sptensor([obj.NumOptParameters, size(obj.ControlPointsRotor)]);
                for iParam = 1:numel(obj.OptimizationParameters)
                    % perturbate parameters
                    Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                    ControlPointsPerturbate = sparse(size(ControlPointsInitRotor, 1), size(ControlPointsInitRotor, 2));
                    optsPerturbate = opts;
                    optsPerturbate.(obj.OptimizationParameters(iParam).Name) = optsPerturbate.(obj.OptimizationParameters(iParam).Name) + Step;
                    [srfPerturbate, ~, ~] = obj.RotorGeometryFunction(optsPerturbate);
                    % See how control points change when parameter changes
                    for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                        ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                        ControlPointsPerturbate(ind_loc, 1) = reshape(srfPerturbate(i).coefs(1, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 2) = reshape(srfPerturbate(i).coefs(2, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 3) = reshape(srfPerturbate(i).coefs(3, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 4) = reshape(srfPerturbate(i).coefs(4, :, :, :), [], 1);
                    end
                    for iCtrPt = 1:numel(ControlPointOffsetsRotor)
                        valsx = ControlPointOffsetsRotor(iCtrPt).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                        valsy = ControlPointOffsetsRotor(iCtrPt).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                        numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                        ControlPointsPerturbate(numbers, 1) = ControlPointsPerturbate(numbers, 1) + valsx;
                        ControlPointsPerturbate(numbers, 2) = ControlPointsPerturbate(numbers, 2) + valsy;
                    end
                    dCrtdP(iParam, :, :) = sptensor(ControlPointsPerturbate - ControlPointsInitRotor)/Step;
                end
                dCstdP = sptensor([obj.NumOptParameters, size(obj.ControlPointsStator)]);
                for iParam = 1:numel(obj.OptimizationParameters)
                    % perturbate parameters
                    Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                    ControlPointsPerturbate = sparse(size(ControlPointsInitStator, 1), size(ControlPointsInitStator, 2));
                    optsPerturbate = opts;
                    optsPerturbate.(obj.OptimizationParameters(iParam).Name) = optsPerturbate.(obj.OptimizationParameters(iParam).Name) + Step;
                    [srfPerturbate, ~, ~] = obj.StatorGeometryFunction(optsPerturbate);
                    % See how control points change when parameter changes
                    for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                        ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                        ControlPointsPerturbate(ind_loc, 1) = reshape(srfPerturbate(i).coefs(1, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 2) = reshape(srfPerturbate(i).coefs(2, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 3) = reshape(srfPerturbate(i).coefs(3, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 4) = reshape(srfPerturbate(i).coefs(4, :, :, :), [], 1);
                    end
                    for iCtrPt = 1:numel(ControlPointOffsetsStator)
                        valsx = ControlPointOffsetsStator(iCtrPt).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                        valsy = ControlPointOffsetsStator(iCtrPt).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                        numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                        ControlPointsPerturbate(numbers, 1) = ControlPointsPerturbate(numbers, 1) + valsx;
                        ControlPointsPerturbate(numbers, 2) = ControlPointsPerturbate(numbers, 2) + valsy;
                    end
                    dCstdP(iParam, :, :) = sptensor(ControlPointsPerturbate - ControlPointsInitStator)/Step;
                end
                % DERIVATIVES stiffness matrix, linear 
                dKrtdClin = sptensor([MultipatchSpaceRotor1.ndof, MultipatchSpaceRotor1.ndof, MultipatchSpaceRotorGeo1.ndof, 2]);
                dKstdClin = sptensor([MultipatchSpaceStator1.ndof, MultipatchSpaceStator1.ndof, MultipatchSpaceStatorGeo1.ndof, 2]);
                for iMat = 1:numel(obj.Materials)
                    if obj.Materials(iMat).Material.getType() == "linear"
                        patchesRt = obj.Materials(iMat).PatchesRotor;
                        patchesSt = obj.Materials(iMat).PatchesStator;
                        nu = obj.Materials(iMat).Material.getNuLinear();
                        dKrtdClin = dKrtdClin + op_D_DC_gradu_gradv_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotor1, spaceRotorEval, ...
                            MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, MultipatchMeshRotor1, mshRotorEval, nu, patchesRt);
                        dKstdClin = dKstdClin + op_D_DC_gradu_gradv_mp_eval(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStator1, spaceStatorEval, ...
                            MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, nu, patchesSt);
                    end
                end

                % derivative force vector wrt control points
                dbrtdC = sptensor([MultipatchSpaceRotor1.ndof, MultipatchSpaceRotorGeo1.ndof, 2]);
                % derivative force vector wrt magnetization of remanent patches
                dbrtdP = sparse(MultipatchSpaceRotor1.ndof, obj.NumOptParameters);
                for iRem = 1:numel(obj.Remanence)
                    remParams = zeros(size(obj.Remanence(iRem).Description));
                    for iRemParam  = 1:numel(obj.Remanence(iRem).Description)
                        remParams(iRemParam) = opts.(obj.Remanence(iRem).Description(iRemParam));
                    end
                    alpha = obj.Remanence(iRem).Relation(remParams);
                    if isempty(alpha)
                        alpha = obj.Remanence(iRem).InitAngle;
                    end
                    patch = obj.Remanence(iRem).Patch;
                    Hc = obj.Remanence(iRem).Material.Br / obj.Remanence(iRem).Material.getMuLinear();
                    % Derivative of control points
                    dbrtdC = dbrtdC + Hc*op_D_DC_gradv_n_bot_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, ...
                        MultipatchMeshRotor1, mshRotorEval, alpha, patch);
                    % Derivative of magnetization parameters, loop if more parameters influence angle
                    dbdAlpha = Hc*op_d_dalpha_gradv_n_bot_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, mshRotorEval, alpha, patch);
                    for iRemParam  = 1:numel(obj.Remanence(iRem).Description)
                        iParam = find(vertcat(obj.OptimizationParameters.Name) == obj.Remanence(iRem).Description(iRemParam));
                        if isempty(iParam)
                            continue
                        end
                        Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                        remParamsPerturbed = remParams;
                        remParamsPerturbed(iRemParam) = remParamsPerturbed(iRemParam) + Step;
                        alpha1 = obj.Remanence(iRem).Relation(remParamsPerturbed);
                        dAlphadP = (alpha1 - alpha)/Step;
%                         dAlphadP = (obj.OptimizationParameters(iParam).MagnetRelations{iMagnetPatch}(x(iParam) + Step) - obj.OptimizationParameters(iParam).MagnetRelations{iMagnetPatch}(x(iParam)))/Step;
                        dbrtdP(:, iParam) = dbrtdP(:, iParam) - dAlphadP*dbdAlpha;
                    end
                end
            end

            % generate linear stiffness matrices
            StiffMatRotor1 = sparse(MultipatchSpaceRotor1.ndof, MultipatchSpaceRotor1.ndof);
            StiffMatStator1 = sparse(MultipatchSpaceStator1.ndof, MultipatchSpaceStator1.ndof);
            for iMat = 1:numel(obj.Materials)
                % here always linear calculation
                patchesRt = obj.Materials(iMat).PatchesRotor;
                patchesSt = obj.Materials(iMat).PatchesStator;
                nu = obj.Materials(iMat).Material.getNuLinear();
                obj.StiffMatsRotor{iMat} = op_gradu_gradv_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotor1, spaceRotorEval, ...
                    MultipatchMeshRotor1, mshRotorEval, nu, patchesRt);
                StiffMatRotor1 = StiffMatRotor1 + obj.StiffMatsRotor{iMat};
                obj.StiffMatsStator{iMat} = op_gradu_gradv_mp_eval(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStator1, spaceStatorEval, ...
                    MultipatchMeshStator1, mshStatorEval, nu, patchesSt);
                StiffMatStator1 = StiffMatStator1 + obj.StiffMatsStator{iMat};
            end
   
            if isfield(opts, 'OPERATING_ANGLE')
                obj.setCurrentOptions(obj.ApplicationCurrent, obj.PolePairs, obj.NumWindings, opts.OPERATING_ANGLE);
            end
%             LENGTH = obj.MachineProperties.Length;
            if isfield(opts, 'LENGTH')
                LENGTH = opts.LENGTH;
            end

            % Iterate over multiple angles and currents
            currents = linspace(0, obj.ApplicationCurrent, 3);
%             currents = obj.ApplicationCurrent;
            angles = obj.OptimizationAngles;
            
%             Torques = cell(numel(currents), numel(angles));
%             dTorques_dP = cell(numel(currents), numel(angles));
%             dTorques_dC = cell(numel(currents), numel(angles));
            Torques = zeros(numel(currents), numel(angles));
            dTorques_dP = zeros(obj.NumOptParameters, numel(currents), numel(angles));
            dTorques_dC = zeros(obj.NumOptControlPointsRotor+obj.NumOptControlPointsStator, numel(currents), numel(angles));


            % init values
            MagneticPotentialRotor1 = obj.MagneticPotentialRotor;
            MagneticPotentialStator1 = obj.MagneticPotentialStator;
            CouplingValues1 = obj.CouplingValues;
            for iCurrent = 1:numel(currents)
                for iAngle = 1:numel(angles)
                    obj.setRotationAngle(angles(iAngle));
                    obj.setCurrent(currents(iCurrent), angles(iAngle));
                    % Nonlinear
                    if obj.OptimizationType == "nonlinear"
                    [MagneticPotentialRotor1, MagneticPotentialStator1, CouplingValues1, A] = ...
                        obj.newtonSolver(MagneticPotentialRotor1, MagneticPotentialStator1, CouplingValues1, ...
                        MultipatchSpaceRotor1, spaceRotorEval, MultipatchMeshRotor1, mshRotorEval, ...
                        MultipatchSpaceStator1, spaceStatorEval, MultipatchMeshStator1, mshStatorEval, false);
                    else
                        [MagneticPotentialRotor1, MagneticPotentialStator1, CouplingValues1, A] = ...
                            obj.linearSolver(StiffMatRotor1, StiffMatStator1);
                    end
                    obj.MagneticPotentialRotor = MagneticPotentialRotor1;
                    obj.MagneticPotentialStator = MagneticPotentialStator1;
                    obj.CouplingValues = CouplingValues1;
%                     obj.plotBResulting();
%                     drawnow();
    
                    T = - MagneticPotentialStator1' * obj.CouplingMatrixStatorInit * obj.RotationMatrixDer * CouplingValues1 * LENGTH * ScaleR^2;
                    Torques(iCurrent, iAngle) = T;
        
                    if nargout == 2
                        dT_du1 = - zeros(size(obj.MagneticPotentialRotor));
                        dT_du2 = - obj.CouplingMatrixStatorInit*obj.RotationMatrixDer*CouplingValues1;
                        dT_du3 = - obj.RotationMatrixDer'*obj.CouplingMatrixStatorInit'*MagneticPotentialStator1;
                        
                        dT_du = [dT_du1(obj.IndependentDOFsRotor); ...
                                  dT_du2(obj.IndependentDOFsStator); ...
                                  dT_du3] * LENGTH * ScaleR^2;
                        % Nonlin
                        dKrtdCnonlin = sptensor([MultipatchSpaceRotor1.ndof, MultipatchSpaceRotor1.ndof, MultipatchSpaceRotorGeo1.ndof, 2]);
                        dKstdCnonlin = sptensor([MultipatchSpaceStator1.ndof, MultipatchSpaceStator1.ndof, MultipatchSpaceStatorGeo1.ndof, 2]);
                        for iMat = 1:numel(obj.Materials)
                            if obj.Materials(iMat).Material.getType() == "nonlinear"
                                patchesRt = obj.Materials(iMat).PatchesRotor;
                                dKrtdCnonlin = dKrtdCnonlin + op_D_DC_gradu_nu_gradv_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotor1, spaceRotorEval, ...
                                    MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, MultipatchMeshRotor1, mshRotorEval, MagneticPotentialRotor1, obj.Materials(iMat).Material, patchesRt);
                                patchesSt = obj.Materials(iMat).PatchesStator;
                                dKstdCnonlin = dKstdCnonlin + op_D_DC_gradu_nu_gradv_mp_eval(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStator1, spaceStatorEval, ...
                                    MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, MagneticPotentialStator1, obj.Materials(iMat).Material, patchesSt);
                            end
                        end
                        dKrtdC = dKrtdClin + dKrtdCnonlin;
                        % Antiperiodic boundary conditions
                        dKrtdC(obj.APdofsRotorLeft', :, :, :) = dKrtdC(obj.APdofsRotorLeft', :, :, :) - dKrtdC(obj.APdofsRotorRight', :, :, :);
                        dKrtdC(obj.APdofsRotorRight', :, :, :) = 0;
    
                        dKstdC = dKstdClin + dKstdCnonlin;
                        % Antiperiodic boundary conditions
                        dKstdC(obj.APdofsStatorLeft', :, :, :) = dKstdC(obj.APdofsStatorLeft', :, :, :) - dKstdC(obj.APdofsStatorRight', :, :, :);
                        dKstdC(obj.APdofsStatorRight', :, :, :) = 0;
        
                        lambda = -A'\dT_du;
    
                        dbstdC = sptensor([MultipatchSpaceStator1.ndof, MultipatchSpaceStatorGeo1.ndof, 2]);
                        dbstdP = sparse(MultipatchSpaceStator1.ndof, obj.NumOptParameters);
                        angleRad = deg2rad(angles(iAngle));
                        % RHS Stator update for angle
                        for iPhase = 1:numel(obj.Phases)
                            switch obj.Phases(iPhase).Type
                                case "U"
                                    Current = currents(iCurrent) * sin (obj.Phase0 + obj.PolePairs * angleRad);
                                case "V"
                                    Current = currents(iCurrent) * sin (obj.Phase0 + obj.PolePairs * angleRad - 2/3*pi);
                                case "W"
                                    Current = currents(iCurrent) * sin (obj.Phase0 + obj.PolePairs * angleRad + 2/3*pi);
                                otherwise
                                    error("Unknown phase type");
                            end
                            dbstdC = dbstdC + obj.Phases(iPhase).dbdC * Current * obj.NumWindings * obj.Phases(iPhase).Slot / obj.Phases(iPhase).Area;
                            % THESE THREE LINES FIX HOW changing
                            % parameters/control points also change the current
                            % density (of lower slots) (product rule)
                            temp = full(obj.Phases(iPhase).Area*obj.Phases(1).dAdC - obj.Phases(1).Area*obj.Phases(iPhase).dAdC);
                            temp = reshape(temp, [1, size(temp)]);
                            dbstdC = dbstdC + ScaleR*obj.Phases(iPhase).Slot*Current/obj.ApplicationCurrent*obj.Jinit/(obj.Phases(iPhase).Area^2)*sptensor(temp.*obj.Phases(iPhase).RHS);
                        
                            for iParam = 1:obj.NumOptParameters
                                if obj.OptimizationParameters(iParam).Name == "ScaleR"
                                    dbstdP(:, iParam) = dbstdP(:, iParam) - Current * obj.NumWindings * obj.Phases(iPhase).Slot / obj.Phases(iPhase).Area*obj.Phases(iPhase).RHS/ScaleR;
                                end
                            end
                        end

                        % Control points
                        dTdCrt = double(ttv(ttv(dKrtdC(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor, :, :), MagneticPotentialRotor1(obj.IndependentDOFsRotor), 2) - dbrtdC(obj.IndependentDOFsRotor, :, :), lambda(1:numel(obj.IndependentDOFsRotor)), 1));     % Tensor toolbox Tensor-vector multiplication
                        dTdCrt_wanted = zeros(obj.NumOptControlPointsRotor, 1);
                        for iCtrPt = 1:obj.NumOptControlPointsRotor
                            numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                            dTdCrt_wanted(iCtrPt) = dTdCrt_wanted(iCtrPt) + ...
                                sum(dTdCrt(numbers, 1).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle) + ...
                                    dTdCrt(numbers, 2).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle));
                        end
                        dTdCst = double(ttv(ttv(dKstdC(obj.IndependentDOFsStator, obj.IndependentDOFsStator, :, :), MagneticPotentialStator1(obj.IndependentDOFsStator), 2) - dbstdC(obj.IndependentDOFsStator, :, :), lambda(numel(obj.IndependentDOFsRotor)+(1:numel(obj.IndependentDOFsStator))), 1));     % Tensor toolbox Tensor-vector multiplication
                        dTdCst_wanted = zeros(obj.NumOptControlPointsStator, 1);
                        for iCtrPt = 1:obj.NumOptControlPointsStator
                            numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                            dTdCst_wanted(iCtrPt) = dTdCst_wanted(iCtrPt) + ...
                                sum(dTdCst(numbers, 1).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle) + ...
                                    dTdCst(numbers, 2).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle));
                        end
                        dTorques_dC(:, iCurrent, iAngle) = [dTdCrt_wanted; dTdCst_wanted];
                        % Geometry parameters
                        dTdP = double(ttt(dCrtdP(:, :, 1:2), sptensor(dTdCrt), [2, 3], [1, 2])) + ...
                            double(ttt(dCstdP(:, :, 1:2), sptensor(dTdCst), [2, 3], [1, 2]));
    
                        % OPERATING_ANGLE
                        for iParam = 1:obj.NumOptParameters
                            if obj.OptimizationParameters(iParam).Name == "OPERATING_ANGLE"
                                dTdP(iParam) = -pi/180*obj.dRHSStator(obj.IndependentDOFsStator)'*lambda(numel(obj.IndependentDOFsRotor)+1:numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator));
                            end
                            if obj.OptimizationParameters(iParam).Name == "LENGTH"
                                dTdP(iParam) = T/LENGTH;
                            end
                            if obj.OptimizationParameters(iParam).Name == "ScaleR"
                                dTdP(iParam) = 2*T/ScaleR;
                            end
                        end
    
                        % add influence of parameters that change magnetization
                        dTdP = dTdP + dbrtdP(obj.IndependentDOFsRotor, :)'*lambda(1:numel(obj.IndependentDOFsRotor));
                        % add influence of parameters that change currents
                        dTdP = dTdP + dbstdP(obj.IndependentDOFsStator, :)'*lambda(numel(obj.IndependentDOFsRotor)+(1:numel(obj.IndependentDOFsStator)));
                        
                        dTorques_dP(:, iCurrent, iAngle) = dTdP;
                    end
                end
            end

            % Objective function evaluation
            PowerConst = obj.Jinit^2/obj.NumWindings/59600000;
            Power = PowerConst*LENGTH*ScaleR^2*obj.Phases(1).Area;
            Cost = 0;
            for iMat = 1:numel(obj.Materials)
                patchesRt = obj.Materials(iMat).PatchesRotor;
                patchesSt = obj.Materials(iMat).PatchesStator;
                rho = obj.Materials(iMat).Material.getRho();
                cost = obj.Materials(iMat).Material.getCost();
                Cost = Cost +  rho * cost * LENGTH * ScaleR^2 * op_Omega_mp_eval(mshRotorEval, patchesRt);
                Cost = Cost +  rho * cost * LENGTH * ScaleR^2 * op_Omega_mp_eval(mshStatorEval, patchesSt);
            end
            
            T_mean = mean(Torques, 2);
            T_sigma = std(Torques, 1, 2); % N instead of N-1

            fopt = obj.MultCost*Cost + obj.MultTsig*sum(T_sigma) + obj.MultPower*Power;
            Mult.MultCost = obj.MultCost;
            Mult.MultTsig = obj.MultTsig;
            Mult.MultPower = obj.MultPower;
            % store values
            obj.OptimizationHistory(end+1).funccount = numel(obj.OptimizationHistory)+1;
            obj.OptimizationHistory(end).fval = fopt;
            obj.OptimizationHistory(end).iteration = -1;
            obj.OptimizationHistory(end).fvalUsed = 0; % will be overwritten by storeOutput function
            obj.OptimizationHistory(end).Tmean = T_mean;
            obj.OptimizationHistory(end).Tsigma = T_sigma;
            obj.OptimizationHistory(end).Cost = Cost;
            obj.OptimizationHistory(end).Currents = currents;
            obj.OptimizationHistory(end).Phases = obj.Phases;
            obj.OptimizationHistory(end).Torques = Torques;
            obj.OptimizationHistory(end).Angles = obj.OptimizationAngles;
            obj.OptimizationHistory(end).Power = Power;
            obj.OptimizationHistory(end).TimeElapsed = toc(obj.TimeStart);
            obj.OptimizationHistory(end).Mult = Mult;
            
            % Postprocessing for mean and std torque derivatives
            if nargout == 2
                dTdP_mean = mean(dTorques_dP, 3);
                dTdC_mean = mean(dTorques_dC, 3);

                % Magnet mass
                gradCostRt = zeros(MultipatchSpaceRotorGeo1.ndof, 2);
                gradCostSt = zeros(MultipatchSpaceStatorGeo1.ndof, 2);
                for iMat = 1:numel(obj.Materials)
                    patchesRt = obj.Materials(iMat).PatchesRotor;
                    patchesSt = obj.Materials(iMat).PatchesStator;
                    rho = obj.Materials(iMat).Material.getRho();
                    cost = obj.Materials(iMat).Material.getCost();
                    gradCostRt = gradCostRt + rho*cost*LENGTH*ScaleR^2*op_D_DC_Omega_mp_eval(MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, MultipatchMeshRotor1, mshRotorEval, patchesRt);
                    gradCostSt = gradCostSt + rho*cost*LENGTH*ScaleR^2*op_D_DC_Omega_mp_eval(MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, patchesSt);
                end
                gradCostP = obj.MultCost*double(ttt(sptensor(dCrtdP(:, :, 1:2)), sptensor(gradCostRt), [2, 3], [1, 2])) + ...
                            obj.MultCost*double(ttt(sptensor(dCstdP(:, :, 1:2)), sptensor(gradCostSt), [2, 3], [1, 2]));

                gradPowerC = PowerConst*LENGTH*ScaleR^2*op_D_DC_Omega_mp_eval(MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, obj.Phases(1).Patches);
                gradPowerP = obj.MultPower* double(ttt(sptensor(dCstdP(:, :, 1:2)), sptensor(gradPowerC), [2, 3], [1, 2]));

                for iParam = 1:obj.NumOptParameters
                    if obj.OptimizationParameters(iParam).Name == "LENGTH"
                        gradCostP(iParam) = gradCostP(iParam) + obj.MultCost*Cost/LENGTH;
                        gradPowerP(iParam) = gradPowerP(iParam) + obj.MultPower*Power/LENGTH;
                    end
                    if obj.OptimizationParameters(iParam).Name == "ScaleR"
                        gradCostP(iParam) = gradCostP(iParam) + 2*obj.MultCost*Cost/ScaleR;
                        gradPowerP(iParam) = gradPowerP(iParam) + 2*obj.MultPower*Power/ScaleR;
                    end
                end

                gradCostCrtwanted = zeros(obj.NumOptControlPointsRotor, 1);
                for iCtrPt = 1:obj.NumOptControlPointsRotor
                    numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                    gradCostCrtwanted(iCtrPt) = gradCostCrtwanted(iCtrPt) + ...
                        sum(gradCostRt(numbers, 1).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle) + ...
                            gradCostRt(numbers, 2).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle));
                end
                gradCostCstwanted = zeros(obj.NumOptControlPointsStator, 1);
                for iCtrPt = 1:obj.NumOptControlPointsStator
                    numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                    gradCostCstwanted(iCtrPt) = gradCostCstwanted(iCtrPt) + ...
                        sum(gradCostSt(numbers, 1).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle) + ...
                            gradCostSt(numbers, 2).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle));
                end
                gradCostC = obj.MultCost*[gradCostCrtwanted; gradCostCstwanted];

%                 gradTsigP = obj.MultTsig/T_sigma*(1/numel(angles)*sum(horzcat(Torques{:}).*horzcat(dTorques_dP{:}), 2) - T_mean*dTdP_mean);
%                 gradTsigC = obj.MultTsig/T_sigma*(1/numel(angles)*sum(horzcat(Torques{:}).*horzcat(dTorques_dC{:}), 2) - T_mean*dTdC_mean);
                gradTsigP = zeros(obj.NumOptParameters, 1);
                gradTsigC = zeros(obj.NumOptControlPointsRotor+obj.NumOptControlPointsStator, 1);
                for iCurrent = 1:numel(currents)
                    gradTsigP = gradTsigP + obj.MultTsig/T_sigma(iCurrent)*(1/numel(angles)*sum(Torques(iCurrent,:).*squeeze(dTorques_dP(:,iCurrent,:)), 2) - T_mean(iCurrent).*dTdP_mean(:, iCurrent));
                    gradTsigC = gradTsigC + obj.MultTsig/T_sigma(iCurrent)*(1/numel(angles)*sum(Torques(iCurrent,:).*squeeze(dTorques_dC(:,iCurrent,:)), 2) - T_mean(iCurrent).*dTdC_mean(:, iCurrent));
                end

                gradfP = gradTsigP + gradCostP + gradPowerP;
                gradfC = gradTsigC + gradCostC;
%                 gradfP = gradfP.*(vertcat(obj.OptimizationParameters.MaxVal) - vertcat(obj.OptimizationParameters.MinVal));
%                 gradfC = gradfC.*([vertcat(obj.OptimizationControlPointsRotor.MaxOffset); vertcat(obj.OptimizationControlPointsStator.MaxOffset)] - [vertcat(obj.OptimizationControlPointsRotor.MinOffset); vertcat(obj.OptimizationControlPointsStator.MinOffset)]);
                gradf = (obj.UpperBounds-obj.LowerBounds).*[gradfP; gradfC];

                obj.LastX = x;
                obj.LastT = mean(Torques(end,:));
                obj.LastdT = (obj.UpperBounds-obj.LowerBounds).*[dTdP_mean(:,end); dTdC_mean(:,end)];
            end
        end

        function [c, gradc] = torqueConstraint(obj, x)
            % Rescale to real values
            x = x.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;
            % If value has already been calculated in objective function:
            if all(x==obj.LastX)
                c = - obj.LastT + obj.MinOptiTorque;
                gradc = - obj.LastdT;
                return;
            end

            opts = obj.Options;
            opts.draw_geometry = false;
            for iParam  = 1:numel (obj.OptimizationParameters)
              opts.(obj.OptimizationParameters(iParam).Name) = x(iParam);
            end
            % Surface and control points of initial geometry
            % Rotor
            [srfInitRotor, ~, ~] = obj.RotorGeometryFunction(opts);
            ControlPointsInitRotor = sparse(size(obj.ControlPointsRotor, 1), size(obj.ControlPointsRotor, 2));
            for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                ControlPointsInitRotor(ind_loc, 1) = reshape(srfInitRotor(i).coefs(1, :, :, :)./srfInitRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitRotor(ind_loc, 2) = reshape(srfInitRotor(i).coefs(2, :, :, :)./srfInitRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitRotor(ind_loc, 3) = reshape(srfInitRotor(i).coefs(3, :, :, :)./srfInitRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitRotor(ind_loc, 4) = reshape(srfInitRotor(i).coefs(4, :, :, :), [], 1);
            end
            % Stator
            [srfInitStator, ~, ~] = obj.StatorGeometryFunction(opts);
            ControlPointsInitStator = sparse(size(obj.ControlPointsStator, 1), size(obj.ControlPointsStator, 2));
            for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                ControlPointsInitStator(ind_loc, 1) = reshape(srfInitStator(i).coefs(1, :, :, :)./srfInitStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitStator(ind_loc, 2) = reshape(srfInitStator(i).coefs(2, :, :, :)./srfInitStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitStator(ind_loc, 3) = reshape(srfInitStator(i).coefs(3, :, :, :)./srfInitStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInitStator(ind_loc, 4) = reshape(srfInitStator(i).coefs(4, :, :, :), [], 1);
            end
            
            % Apply control point offsets
            ControlPointOffsetsRotor = x(obj.NumOptParameters+(1:obj.NumOptControlPointsRotor));
            for iCtrPt = 1:numel(ControlPointOffsetsRotor)
                valsx = ControlPointOffsetsRotor(iCtrPt).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                valsy = ControlPointOffsetsRotor(iCtrPt).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                ControlPointsInitRotor(numbers, 1) = ControlPointsInitRotor(numbers, 1) + valsx;
                ControlPointsInitRotor(numbers, 2) = ControlPointsInitRotor(numbers, 2) + valsy;
            end
            % Stator one
            ControlPointOffsetsStator = x(obj.NumOptParameters+obj.NumOptControlPointsRotor+(1:obj.NumOptControlPointsStator));
            for iCtrPt = 1:numel(ControlPointOffsetsStator)
                valsx = ControlPointOffsetsStator(iCtrPt).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                valsy = ControlPointOffsetsStator(iCtrPt).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                ControlPointsInitStator(numbers, 1) = ControlPointsInitStator(numbers, 1) + valsx;
                ControlPointsInitStator(numbers, 2) = ControlPointsInitStator(numbers, 2) + valsy;
            end
            % Update initial surface
            % Rotor
            for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                srfInitRotor(i).coefs(1, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 1).*ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
                srfInitRotor(i).coefs(2, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 2).*ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
                srfInitRotor(i).coefs(3, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 3).*ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
                srfInitRotor(i).coefs(4, :, :, :) = reshape(ControlPointsInitRotor(ind_loc, 4), srfInitRotor(i).number);
            end
            % Stator
            for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                srfInitStator(i).coefs(1, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 1).*ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
                srfInitStator(i).coefs(2, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 2).*ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
                srfInitStator(i).coefs(3, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 3).*ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
                srfInitStator(i).coefs(4, :, :, :) = reshape(ControlPointsInitStator(ind_loc, 4), srfInitStator(i).number);
            end
            
            % Load modified geometries: Rotor
            [GeometryRotor1, ~, ~, ~, ~] = mp_geo_load(srfInitRotor);
            meshes = cell (1, obj.NumberPatchesRotor);
            spaces  = cell (1, obj.NumberPatchesRotor);
            spaces_geo = cell (1, obj.NumberPatchesRotor);
            % create space partition and form functions
            for iptc = 1:obj.NumberPatchesRotor
              [knots, zeta] = kntrefine (GeometryRotor1(iptc).nurbs.knots, obj.SubdivisionsPatchesRotor-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
              rule      = msh_gauss_nodes (obj.DegreeQuadrature);
              [qn, qw]  = msh_set_quad_nodes (zeta, rule);
              meshes{iptc} = msh_cartesian (zeta, qn, qw, GeometryRotor1(iptc));
              spaces{iptc}  = sp_bspline (knots, obj.FormFunctionDegree, meshes{iptc});
              spaces_geo{iptc} = sp_nurbs(GeometryRotor1(iptc).nurbs, meshes{iptc});
            end
            MultipatchMeshRotor1 = msh_multipatch(meshes, obj.BoundariesRotor);
            MultipatchSpaceRotor1  = sp_multipatch(spaces, MultipatchMeshRotor1, obj.InterfacesRotor, obj.BoundaryInterfacesRotor);
            MultipatchSpaceRotorGeo1 = sp_multipatch(spaces_geo, MultipatchMeshRotor1, obj.InterfacesRotor, obj.BoundaryInterfacesRotor);
            % Precalculate Spaces and Meshes
            mshRotorEval = cell (1, obj.NumberPatchesRotor);
            spaceRotorEval  = cell (1, obj.NumberPatchesRotor);
            spaceRotorEvalGeo = cell (1, obj.NumberPatchesRotor);
            for iPatch = 1:obj.NumberPatchesRotor
                mshRotorEval{iPatch} = msh_precompute(MultipatchMeshRotor1.msh_patch{iPatch});
                spaceRotorEval{iPatch}  = sp_precompute(MultipatchSpaceRotor1.sp_patch{iPatch}, mshRotorEval{iPatch}, 'gradient', true);
                spaceRotorEvalGeo{iPatch}  = sp_precompute_param(MultipatchSpaceRotorGeo1.sp_patch{iPatch}, mshRotorEval{iPatch}, 'value', true, 'gradient', true);
            end

            [GeometryStator1, ~, ~, ~, ~] = mp_geo_load(srfInitStator);
            meshes = cell (1, obj.NumberPatchesStator);
            spaces  = cell (1, obj.NumberPatchesStator);
            spaces_geo = cell (1, obj.NumberPatchesStator);
            % create space partition and form functions
            for iptc = 1:obj.NumberPatchesStator
              [knots, zeta] = kntrefine (GeometryStator1(iptc).nurbs.knots, obj.SubdivisionsPatchesStator-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
              rule      = msh_gauss_nodes (obj.DegreeQuadrature);
              [qn, qw]  = msh_set_quad_nodes (zeta, rule);
              meshes{iptc} = msh_cartesian (zeta, qn, qw, GeometryStator1(iptc));
              spaces{iptc}  = sp_bspline (knots, obj.FormFunctionDegree, meshes{iptc});
              spaces_geo{iptc} = sp_nurbs(GeometryStator1(iptc).nurbs, meshes{iptc});
            end
            MultipatchMeshStator1 = msh_multipatch(meshes, obj.BoundariesStator);
            MultipatchSpaceStator1  = sp_multipatch(spaces, MultipatchMeshStator1, obj.InterfacesStator, obj.BoundaryInterfacesStator);
            MultipatchSpaceStatorGeo1 = sp_multipatch(spaces_geo, MultipatchMeshStator1, obj.InterfacesStator, obj.BoundaryInterfacesStator);
            % Precalculate Spaces and Meshes
            mshStatorEval = cell (1, obj.NumberPatchesStator);
            spaceStatorEval  = cell (1, obj.NumberPatchesStator);
            spaceStatorEvalGeo = cell (1, obj.NumberPatchesStator);
            for iPatch = 1:obj.NumberPatchesStator
                mshStatorEval{iPatch} = msh_precompute(MultipatchMeshStator1.msh_patch{iPatch});
                spaceStatorEval{iPatch}  = sp_precompute(MultipatchSpaceStator1.sp_patch{iPatch}, mshStatorEval{iPatch}, 'gradient', true);
                spaceStatorEvalGeo{iPatch}  = sp_precompute_param(MultipatchSpaceStatorGeo1.sp_patch{iPatch}, mshStatorEval{iPatch}, 'value', true, 'gradient', true);
            end

            ScaleR = 1;
            if isfield(opts, 'ScaleR')
                ScaleR = opts.ScaleR;
            end

            % Build RHS for rotor with updated magnetization depending on parameters that change magnetization direction
            RHS1 = sparse(MultipatchSpaceRotor1.ndof, 1);
            for iRem = 1:numel(obj.Remanence)
                remParams = zeros(size(obj.Remanence(iRem).Description));
                for iParam  = 1:numel(obj.Remanence(iRem).Description)
                    remParams(iParam) = opts.(obj.Remanence(iRem).Description(iParam));
                end
                alpha = obj.Remanence(iRem).Relation(remParams);
                if isempty(alpha)
                    alpha = obj.Remanence(iRem).InitAngle;
                end
                patch = obj.Remanence(iRem).Patch;
                Hc = obj.Remanence(iRem).Material.Br / obj.Remanence(iRem).Material.getMuLinear();
                RHS1 = RHS1 + Hc*op_gradv_n_bot_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, mshRotorEval, alpha, patch);
            end
            obj.RHSRotor = RHS1;

            % Reconstruct RHS of stator for updated geometry
            obj.optimizationGeneratePhases(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval,ScaleR);

            if nargout == 2
                % Numerical derivatives dC/dP of parameter iOpt
                dCrtdP = sptensor([obj.NumOptParameters, size(obj.ControlPointsRotor)]);
                for iParam = 1:numel(obj.OptimizationParameters)
                    % perturbate parameters
                    Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                    ControlPointsPerturbate = sparse(size(ControlPointsInitRotor, 1), size(ControlPointsInitRotor, 2));
                    optsPerturbate = opts;
                    optsPerturbate.(obj.OptimizationParameters(iParam).Name) = optsPerturbate.(obj.OptimizationParameters(iParam).Name) + Step;
                    [srfPerturbate, ~, ~] = obj.RotorGeometryFunction(optsPerturbate);
                    % See how control points change when parameter changes
                    for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                        ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                        ControlPointsPerturbate(ind_loc, 1) = reshape(srfPerturbate(i).coefs(1, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 2) = reshape(srfPerturbate(i).coefs(2, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 3) = reshape(srfPerturbate(i).coefs(3, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 4) = reshape(srfPerturbate(i).coefs(4, :, :, :), [], 1);
                    end
                    for iCtrPt = 1:numel(ControlPointOffsetsRotor)
                        valsx = ControlPointOffsetsRotor(iCtrPt).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                        valsy = ControlPointOffsetsRotor(iCtrPt).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                        numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                        ControlPointsPerturbate(numbers, 1) = ControlPointsPerturbate(numbers, 1) + valsx;
                        ControlPointsPerturbate(numbers, 2) = ControlPointsPerturbate(numbers, 2) + valsy;
                    end
                    dCrtdP(iParam, :, :) = sptensor(ControlPointsPerturbate - ControlPointsInitRotor)/Step;
                end
                dCstdP = sptensor([obj.NumOptParameters, size(obj.ControlPointsStator)]);
                for iParam = 1:numel(obj.OptimizationParameters)
                    % perturbate parameters
                    Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                    ControlPointsPerturbate = sparse(size(ControlPointsInitStator, 1), size(ControlPointsInitStator, 2));
                    optsPerturbate = opts;
                    optsPerturbate.(obj.OptimizationParameters(iParam).Name) = optsPerturbate.(obj.OptimizationParameters(iParam).Name) + Step;
                    [srfPerturbate, ~, ~] = obj.StatorGeometryFunction(optsPerturbate);
                    % See how control points change when parameter changes
                    for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                        ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                        ControlPointsPerturbate(ind_loc, 1) = reshape(srfPerturbate(i).coefs(1, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 2) = reshape(srfPerturbate(i).coefs(2, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 3) = reshape(srfPerturbate(i).coefs(3, :, :, :)./srfPerturbate(i).coefs(4, :, :, :), [], 1);
                        ControlPointsPerturbate(ind_loc, 4) = reshape(srfPerturbate(i).coefs(4, :, :, :), [], 1);
                    end
                    for iCtrPt = 1:numel(ControlPointOffsetsStator)
                        valsx = ControlPointOffsetsStator(iCtrPt).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                        valsy = ControlPointOffsetsStator(iCtrPt).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                        numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                        ControlPointsPerturbate(numbers, 1) = ControlPointsPerturbate(numbers, 1) + valsx;
                        ControlPointsPerturbate(numbers, 2) = ControlPointsPerturbate(numbers, 2) + valsy;
                    end
                    dCstdP(iParam, :, :) = sptensor(ControlPointsPerturbate - ControlPointsInitStator)/Step;
                end
                % DERIVATIVES stiffness matrix, linear 
                dKrtdClin = sptensor([MultipatchSpaceRotor1.ndof, MultipatchSpaceRotor1.ndof, MultipatchSpaceRotorGeo1.ndof, 2]);
                dKstdClin = sptensor([MultipatchSpaceStator1.ndof, MultipatchSpaceStator1.ndof, MultipatchSpaceStatorGeo1.ndof, 2]);
                for iMat = 1:numel(obj.Materials)
                    if obj.Materials(iMat).Material.getType() == "linear"
                        patchesRt = obj.Materials(iMat).PatchesRotor;
                        patchesSt = obj.Materials(iMat).PatchesStator;
                        nu = obj.Materials(iMat).Material.getNuLinear();
                        dKrtdClin = dKrtdClin + op_D_DC_gradu_gradv_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotor1, spaceRotorEval, ...
                            MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, MultipatchMeshRotor1, mshRotorEval, nu, patchesRt);
                        dKstdClin = dKstdClin + op_D_DC_gradu_gradv_mp_eval(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStator1, spaceStatorEval, ...
                            MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, nu, patchesSt);
                    end
                end

                % derivative force vector wrt control points
                dbrtdC = sptensor([MultipatchSpaceRotor1.ndof, MultipatchSpaceRotorGeo1.ndof, 2]);
                % derivative force vector wrt magnetization of remanent patches
                dbdP = sparse(MultipatchSpaceRotor1.ndof, obj.NumOptParameters);
                for iRem = 1:numel(obj.Remanence)
                    remParams = zeros(size(obj.Remanence(iRem).Description));
                    for iRemParam  = 1:numel(obj.Remanence(iRem).Description)
                        remParams(iRemParam) = opts.(obj.Remanence(iRem).Description(iRemParam));
                    end
                    alpha = obj.Remanence(iRem).Relation(remParams);
                    if isempty(alpha)
                        alpha = obj.Remanence(iRem).InitAngle;
                    end
                    patch = obj.Remanence(iRem).Patch;
                    Hc = obj.Remanence(iRem).Material.Br / obj.Remanence(iRem).Material.getMuLinear();
                    % Derivative of control points
                    dbrtdC = dbrtdC + Hc*op_D_DC_gradv_n_bot_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, ...
                        MultipatchMeshRotor1, mshRotorEval, alpha, patch);
                    % Derivative of magnetization parameters, loop if more parameters influence angle
                    dbdAlpha = Hc*op_d_dalpha_gradv_n_bot_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, mshRotorEval, alpha, patch);
                    for iRemParam  = 1:numel(obj.Remanence(iRem).Description)
                        iParam = find(vertcat(obj.OptimizationParameters.Name) == obj.Remanence(iRem).Description(iRemParam));
                        if isempty(iParam)
                            continue
                        end
                        Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                        remParamsPerturbed = remParams;
                        remParamsPerturbed(iRemParam) = remParamsPerturbed(iRemParam) + Step;
                        alpha1 = obj.Remanence(iRem).Relation(remParamsPerturbed);
                        dAlphadP = (alpha1 - alpha)/Step;
%                         dAlphadP = (obj.OptimizationParameters(iParam).MagnetRelations{iMagnetPatch}(x(iParam) + Step) - obj.OptimizationParameters(iParam).MagnetRelations{iMagnetPatch}(x(iParam)))/Step;
                        dbdP(:, iParam) = dbdP(:, iParam) + dAlphadP*dbdAlpha;
                    end
                end
            end

            % generate linear stiffness matrices
            StiffMatRotor1 = sparse(MultipatchSpaceRotor1.ndof, MultipatchSpaceRotor1.ndof);
            StiffMatStator1 = sparse(MultipatchSpaceStator1.ndof, MultipatchSpaceStator1.ndof);
            for iMat = 1:numel(obj.Materials)
                % here always linear calculation
                patchesRt = obj.Materials(iMat).PatchesRotor;
                patchesSt = obj.Materials(iMat).PatchesStator;
                nu = obj.Materials(iMat).Material.getNuLinear();
                obj.StiffMatsRotor{iMat} = op_gradu_gradv_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotor1, spaceRotorEval, ...
                    MultipatchMeshRotor1, mshRotorEval, nu, patchesRt);
                StiffMatRotor1 = StiffMatRotor1 + obj.StiffMatsRotor{iMat};
                obj.StiffMatsStator{iMat} = op_gradu_gradv_mp_eval(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStator1, spaceStatorEval, ...
                    MultipatchMeshStator1, mshStatorEval, nu, patchesSt);
                StiffMatStator1 = StiffMatStator1 + obj.StiffMatsStator{iMat};
            end
   
            if isfield(opts, 'OPERATING_ANGLE')
                obj.setCurrentOptions(obj.ApplicationCurrent, obj.PolePairs, obj.NumWindings, opts.OPERATING_ANGLE);
            end
%             LENGTH = obj.MachineProperties.Length;
            if isfield(opts, 'LENGTH')
                LENGTH = opts.LENGTH;
            end

            % Iterate over multiple angles and currents
            currents = linspace(0, obj.ApplicationCurrent, 3);
%             currents = obj.ApplicationCurrent;
            angles = obj.OptimizationAngles;
            
            Torques = zeros(numel(currents), numel(angles));
            dTorques_dP = zeros(obj.NumOptParameters, numel(currents), numel(angles));
            dTorques_dC = zeros(obj.NumOptControlPointsRotor+obj.NumOptControlPointsStator, numel(currents), numel(angles));


            % init values
            MagneticPotentialRotor1 = obj.MagneticPotentialRotor;
            MagneticPotentialStator1 = obj.MagneticPotentialStator;
            CouplingValues1 = obj.CouplingValues;
            for iCurrent = numel(currents)
                for iAngle = 1:numel(angles)
                    obj.setRotationAngle(angles(iAngle));
                    obj.setCurrent(currents(iCurrent), angles(iAngle));
                    % Nonlinear
                    if obj.OptimizationType == "nonlinear"
                    [MagneticPotentialRotor1, MagneticPotentialStator1, CouplingValues1, A] = ...
                        obj.newtonSolver(MagneticPotentialRotor1, MagneticPotentialStator1, CouplingValues1, ...
                        MultipatchSpaceRotor1, spaceRotorEval, MultipatchMeshRotor1, mshRotorEval, ...
                        MultipatchSpaceStator1, spaceStatorEval, MultipatchMeshStator1, mshStatorEval, false);
                    else
                        [MagneticPotentialRotor1, MagneticPotentialStator1, CouplingValues1, A] = ...
                            obj.linearSolver(StiffMatRotor1, StiffMatStator1);
                    end
                    obj.MagneticPotentialRotor = MagneticPotentialRotor1;
                    obj.MagneticPotentialStator = MagneticPotentialStator1;
                    obj.CouplingValues = CouplingValues1;
%                     obj.plotBResulting();
%                     drawnow();
    
                    T = - MagneticPotentialStator1' * obj.CouplingMatrixStatorInit * obj.RotationMatrixDer * CouplingValues1 * LENGTH * ScaleR^2;
                    Torques(iCurrent, iAngle) = T;
        
                    if nargout == 2
                        dT_du1 = - zeros(size(obj.MagneticPotentialRotor));
                        dT_du2 = - obj.CouplingMatrixStatorInit*obj.RotationMatrixDer*CouplingValues1;
                        dT_du3 = - obj.RotationMatrixDer'*obj.CouplingMatrixStatorInit'*MagneticPotentialStator1;
                        
                        dT_du = [dT_du1(obj.IndependentDOFsRotor); ...
                                  dT_du2(obj.IndependentDOFsStator); ...
                                  dT_du3] * LENGTH * ScaleR^2;
                        % Nonlin
                        dKrtdCnonlin = sptensor([MultipatchSpaceRotor1.ndof, MultipatchSpaceRotor1.ndof, MultipatchSpaceRotorGeo1.ndof, 2]);
                        dKstdCnonlin = sptensor([MultipatchSpaceStator1.ndof, MultipatchSpaceStator1.ndof, MultipatchSpaceStatorGeo1.ndof, 2]);
                        for iMat = 1:numel(obj.Materials)
                            if obj.Materials(iMat).Material.getType() == "nonlinear"
                                patchesRt = obj.Materials(iMat).PatchesRotor;
                                dKrtdCnonlin = dKrtdCnonlin + op_D_DC_gradu_nu_gradv_mp_eval(MultipatchSpaceRotor1, spaceRotorEval, MultipatchSpaceRotor1, spaceRotorEval, ...
                                    MultipatchSpaceRotorGeo1, spaceRotorEvalGeo, MultipatchMeshRotor1, mshRotorEval, MagneticPotentialRotor1, obj.Materials(iMat).Material, patchesRt);
                                patchesSt = obj.Materials(iMat).PatchesStator;
                                dKstdCnonlin = dKstdCnonlin + op_D_DC_gradu_nu_gradv_mp_eval(MultipatchSpaceStator1, spaceStatorEval, MultipatchSpaceStator1, spaceStatorEval, ...
                                    MultipatchSpaceStatorGeo1, spaceStatorEvalGeo, MultipatchMeshStator1, mshStatorEval, MagneticPotentialStator1, obj.Materials(iMat).Material, patchesSt);
                            end
                        end
                        dKrtdC = dKrtdClin + dKrtdCnonlin;
                        % Antiperiodic boundary conditions
                        dKrtdC(obj.APdofsRotorLeft', :, :, :) = dKrtdC(obj.APdofsRotorLeft', :, :, :) - dKrtdC(obj.APdofsRotorRight', :, :, :);
                        dKrtdC(obj.APdofsRotorRight', :, :, :) = 0;
    
                        dKstdC = dKstdClin + dKstdCnonlin;
                        % Antiperiodic boundary conditions
                        dKstdC(obj.APdofsStatorLeft', :, :, :) = dKstdC(obj.APdofsStatorLeft', :, :, :) - dKstdC(obj.APdofsStatorRight', :, :, :);
                        dKstdC(obj.APdofsStatorRight', :, :, :) = 0;
        
                        lambda = -A'\dT_du;
    
                        dbstdC = sptensor([MultipatchSpaceStator1.ndof, MultipatchSpaceStatorGeo1.ndof, 2]);
                        dbstdP = sparse(MultipatchSpaceStator1.ndof, obj.NumOptParameters);
                        angleRad = deg2rad(angles(iAngle));
                        % RHS Stator update for angle
                        for iPhase = 1:numel(obj.Phases)
                            switch obj.Phases(iPhase).Type
                                case "U"
                                    Current = currents(iCurrent) * sin (obj.Phase0 + obj.PolePairs * angleRad);
                                case "V"
                                    Current = currents(iCurrent) * sin (obj.Phase0 + obj.PolePairs * angleRad - 2/3*pi);
                                case "W"
                                    Current = currents(iCurrent) * sin (obj.Phase0 + obj.PolePairs * angleRad + 2/3*pi);
                                otherwise
                                    error("Unknown phase type");
                            end
                            dbstdC = dbstdC + obj.Phases(iPhase).dbdC * Current * obj.NumWindings * obj.Phases(iPhase).Slot / obj.Phases(iPhase).Area;
                            % THESE THREE LINES FIX HOW changing
                            % parameters/control points also change the current
                            % density (of lower slots) (product rule)
                            temp = full(obj.Phases(iPhase).Area*obj.Phases(1).dAdC - obj.Phases(1).Area*obj.Phases(iPhase).dAdC);
                            temp = reshape(temp, [1, size(temp)]);
                            dbstdC = dbstdC + ScaleR*obj.Phases(iPhase).Slot*Current/obj.ApplicationCurrent*obj.Jinit/(obj.Phases(iPhase).Area^2)*sptensor(temp.*obj.Phases(iPhase).RHS);
                        
                            for iParam = 1:obj.NumOptParameters
                                if obj.OptimizationParameters(iParam).Name == "ScaleR"
                                    dbstdP(:, iParam) = dbstdP(:, iParam) - Current * obj.NumWindings * obj.Phases(iPhase).Slot / obj.Phases(iPhase).Area*obj.Phases(iPhase).RHS/ScaleR;
                                end
                            end
                        end
            
                        % Control points
                        dTdCrt = double(ttv(ttv(dKrtdC(obj.IndependentDOFsRotor, obj.IndependentDOFsRotor, :, :), MagneticPotentialRotor1(obj.IndependentDOFsRotor), 2) - dbrtdC(obj.IndependentDOFsRotor, :, :), lambda(1:numel(obj.IndependentDOFsRotor)), 1));     % Tensor toolbox Tensor-vector multiplication
                        dTdCrt_wanted = zeros(obj.NumOptControlPointsRotor, 1);
                        for iCtrPt = 1:obj.NumOptControlPointsRotor
                            numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                            dTdCrt_wanted(iCtrPt) = dTdCrt_wanted(iCtrPt) + ...
                                sum(dTdCrt(numbers, 1).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle) + ...
                                    dTdCrt(numbers, 2).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle));
                        end
                        dTdCst = double(ttv(ttv(dKstdC(obj.IndependentDOFsStator, obj.IndependentDOFsStator, :, :), MagneticPotentialStator1(obj.IndependentDOFsStator), 2) - dbstdC(obj.IndependentDOFsStator, :, :), lambda(numel(obj.IndependentDOFsRotor)+(1:numel(obj.IndependentDOFsStator))), 1));     % Tensor toolbox Tensor-vector multiplication
                        dTdCst_wanted = zeros(obj.NumOptControlPointsStator, 1);
                        for iCtrPt = 1:obj.NumOptControlPointsStator
                            numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                            dTdCst_wanted(iCtrPt) = dTdCst_wanted(iCtrPt) + ...
                                sum(dTdCst(numbers, 1).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle) + ...
                                    dTdCst(numbers, 2).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle));
                        end
                        dTorques_dC(:, iCurrent, iAngle) = [dTdCrt_wanted; dTdCst_wanted];
                        % Geometry parameters
                        dTdP = double(ttt(dCrtdP(:, :, 1:2), sptensor(dTdCrt), [2, 3], [1, 2])) + ...
                            double(ttt(dCstdP(:, :, 1:2), sptensor(dTdCst), [2, 3], [1, 2]));
    
                        % OPERATING_ANGLE
                        for iParam = 1:obj.NumOptParameters
                            if obj.OptimizationParameters(iParam).Name == "OPERATING_ANGLE"
                                dTdP(iParam) = -pi/180*obj.dRHSStator(obj.IndependentDOFsStator)'*lambda(numel(obj.IndependentDOFsRotor)+1:numel(obj.IndependentDOFsRotor)+numel(obj.IndependentDOFsStator));
                            end
                            if obj.OptimizationParameters(iParam).Name == "LENGTH"
                                dTdP(iParam) = T/LENGTH;
                            end
                            if obj.OptimizationParameters(iParam).Name == "ScaleR"
                                dTdP(iParam) = 2*T/ScaleR;
                            end
                        end
    
                        % add influence of parameters that change magnetization
                        dTdP = dTdP + dbdP(obj.IndependentDOFsRotor, :)'*lambda(1:numel(obj.IndependentDOFsRotor));
                        % add influence of parameters that change currents
                        dTdP = dTdP + dbstdP(obj.IndependentDOFsStator, :)'*lambda(numel(obj.IndependentDOFsRotor)+(1:numel(obj.IndependentDOFsStator)));
                        
                        dTorques_dP(:, iCurrent, iAngle) = dTdP;
                    end
                end
            end
            
            T_mean = mean(Torques, 2);
            T_sigma = std(Torques, 1, 2); % N instead of N-1
            
            % Postprocessing for mean and std torque derivatives

            dTdP_mean = mean(dTorques_dP, 3);
            dTdC_mean = mean(dTorques_dC, 3);

            c = - mean(Torques(end,:)) + obj.MinOptiTorque;
            gradc = -  (obj.UpperBounds-obj.LowerBounds).*[dTdP_mean(:,end); dTdC_mean(:,end)];
%             obj.LastX = x;
%             obj.LastT = mean(Torques(end,:));
%             obj.LastdT = (obj.UpperBounds-obj.LowerBounds).*[dTdP_mean(:,end); dTdC_mean(:,end)];

            
        end
        
        function x_opt = optimize(obj)
            obj.OptimizationType = "linear";
            for iMat = 1:numel(obj.Materials)
                if obj.Materials(iMat).Material.getType() == "nonlinear"
                    obj.OptimizationType = "nonlinear";
                end
            end
            obj.OptimizationHistory=struct('iteration', {}, 'funccount', {}, 'fval', {}, 'fvalUsed', {}, 'Tmean', {}, 'Tsigma', {}, 'Cost', {}, 'Currents',  {}, 'Angles', {}, 'Torques', {}, 'Parameters', {}, 'TimeElapsed', {}, 'Power', {}, 'Phases', {}, 'Mult', {});
            obj.StartOptTime = string(datetime('now', 'Format', "yyyy-MM-dd-HH-mm"));
            obj.SaveDirectory = "plots/"+obj.StartOptTime;
            mkdir(obj.SaveDirectory);

            x0Params = (vertcat(obj.OptimizationParameters.InitVal) - vertcat(obj.OptimizationParameters.MinVal))./(vertcat(obj.OptimizationParameters.MaxVal) - vertcat(obj.OptimizationParameters.MinVal));
            x0ControlPointsRt = (vertcat(obj.OptimizationControlPointsRotor.InitOffset) - vertcat(obj.OptimizationControlPointsRotor.MinOffset))./(vertcat(obj.OptimizationControlPointsRotor.MaxOffset) - vertcat(obj.OptimizationControlPointsRotor.MinOffset));
            x0ControlPointsSt = (vertcat(obj.OptimizationControlPointsStator.InitOffset) - vertcat(obj.OptimizationControlPointsStator.MinOffset))./(vertcat(obj.OptimizationControlPointsStator.MaxOffset) - vertcat(obj.OptimizationControlPointsStator.MinOffset));

            problem.options = optimoptions(@fmincon, ... % 'UseParallel', true, ...
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                'Display', 'iter-detailed', ...
                'Algorithm', obj.OptimizationMethod, ...
                'MaxIterations', 200, ...
                'OutputFcn', @(a, b, c) obj.storeOutput(a, b, c), ...
                'PlotFcn', @(a, b, c) obj.plotFunction(a, b, c), ...
                'StepTolerance', 1e-8, ... 
                'OptimalityTolerance', 1e-6, ...
                'ConstraintTolerance', 1e-6, ...
                'ScaleProblem', false, ... % scale manually all values (nonlinear would be not scaled otherwise)
                'HonorBounds', false);   % do not reset design variables that are at the bounds
            problem.objective = @(x) obj.optimizationFunctionParameterShape(x);
            problem.nonlcon = @(x) obj.optimizationConstraints(x);
            problem.x0 = [x0Params; x0ControlPointsRt; x0ControlPointsSt];
            problem.Aeq = [];
            problem.beq = [];
            problem.ub = ones(size(vertcat(obj.OptimizationParameters.MaxVal, obj.OptimizationControlPointsRotor.MaxOffset, obj.OptimizationControlPointsStator.MaxOffset)));
            problem.lb = zeros(size(vertcat(obj.OptimizationParameters.MinVal, obj.OptimizationControlPointsRotor.MinOffset, obj.OptimizationControlPointsStator.MinOffset)));
            problem.solver = 'fmincon';

            obj.TimeStart = tic();

            x_opt = fmincon(problem);

            x_opt = x_opt.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;
        end
        
        %% Optimization plot function
        function stop = plotFunction(obj, x, optimValues, state)
            x = x.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            stop = false;
            opts = obj.Options;
            opts.draw_geometry = false;
            for iParam  = 1:numel (obj.OptimizationParameters)
              opts.(obj.OptimizationParameters(iParam).Name) = x(iParam);
            end

            [srfRotor, ~, ~] = obj.RotorGeometryFunction(opts);
            % Extract and update control points and surface
            ControlPointsInit = obj.ControlPointsRotor;
            for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                ControlPointsInit(ind_loc, 1) = reshape(srfRotor(i).coefs(1, :, :, :)./srfRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInit(ind_loc, 2) = reshape(srfRotor(i).coefs(2, :, :, :)./srfRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInit(ind_loc, 3) = reshape(srfRotor(i).coefs(3, :, :, :)./srfRotor(i).coefs(4, :, :, :), [], 1);
                ControlPointsInit(ind_loc, 4) = reshape(srfRotor(i).coefs(4, :, :, :), [], 1);
            end
            % Control point offsets
            vals = x(obj.NumOptParameters+(1:obj.NumOptControlPointsRotor));
            for iCtrPt = 1:numel(vals)
                valsx = vals(iCtrPt).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                valsy = vals(iCtrPt).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                ControlPointsInit(numbers, 1) = ControlPointsInit(numbers, 1) + valsx;
                ControlPointsInit(numbers, 2) = ControlPointsInit(numbers, 2) + valsy;
            end
            % update initial surface
            for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                srfRotor(i).coefs(1, :, :, :) = reshape(ControlPointsInit(ind_loc, 1).*ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                srfRotor(i).coefs(2, :, :, :) = reshape(ControlPointsInit(ind_loc, 2).*ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                srfRotor(i).coefs(3, :, :, :) = reshape(ControlPointsInit(ind_loc, 3).*ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                srfRotor(i).coefs(4, :, :, :) = reshape(ControlPointsInit(ind_loc, 4), srfRotor(i).number);
            end
            % update plot
            PlotPointsRotor = obj.CalculateGeoPlotPointsSrf(srfRotor);

            [srfStator, ~, ~] = obj.StatorGeometryFunction(opts);
            % Extract and update control points and surface
            ControlPointsInit = obj.ControlPointsStator;
            for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                ControlPointsInit(ind_loc, 1) = reshape(srfStator(i).coefs(1, :, :, :)./srfStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInit(ind_loc, 2) = reshape(srfStator(i).coefs(2, :, :, :)./srfStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInit(ind_loc, 3) = reshape(srfStator(i).coefs(3, :, :, :)./srfStator(i).coefs(4, :, :, :), [], 1);
                ControlPointsInit(ind_loc, 4) = reshape(srfStator(i).coefs(4, :, :, :), [], 1);
            end
            % Control point offsets
            vals = x(obj.NumOptParameters+obj.NumOptControlPointsRotor+(1:obj.NumOptControlPointsStator));
            for iCtrPt = 1:numel(vals)
                valsx = vals(iCtrPt).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                valsy = vals(iCtrPt).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                ControlPointsInit(numbers, 1) = ControlPointsInit(numbers, 1) + valsx;
                ControlPointsInit(numbers, 2) = ControlPointsInit(numbers, 2) + valsy;
            end
            % update initial surface
            for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                srfStator(i).coefs(1, :, :, :) = reshape(ControlPointsInit(ind_loc, 1).*ControlPointsInit(ind_loc, 4), srfStator(i).number);
                srfStator(i).coefs(2, :, :, :) = reshape(ControlPointsInit(ind_loc, 2).*ControlPointsInit(ind_loc, 4), srfStator(i).number);
                srfStator(i).coefs(3, :, :, :) = reshape(ControlPointsInit(ind_loc, 3).*ControlPointsInit(ind_loc, 4), srfStator(i).number);
                srfStator(i).coefs(4, :, :, :) = reshape(ControlPointsInit(ind_loc, 4), srfStator(i).number);
            end
            % update plot
            PlotPointsStator = obj.CalculateGeoPlotPointsSrf(srfStator);
            
            switch state
                case 'init'
                    obj.OptimizationFigure = gca;
                    hold on
                    axis equal
                    % Rotor
                    for iMat = 1:numel(obj.Materials)
                        for patch = obj.Materials(iMat).PatchesRotor
                            obj.srfplotRotor{patch}=surf(PlotPointsRotor{patch}{1, 1}, PlotPointsRotor{patch}{2, 1}, zeros(size(obj.GeoPlotPointsRotor{patch}{1, 1})), ...
                                "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.PlotFaceAlpha);
                        end
                    end
                    for iBnd = 1:numel(obj.BoundariesRotor)
                        patch = obj.BoundariesRotor(iBnd).patches;
                        side = obj.BoundariesRotor(iBnd).faces;
                        obj.lplotRotor{iBnd} = plot(PlotPointsRotor{patch}{1, side+1}, PlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                    end
                    for iInt = 1:numel(obj.InterfacesRotor)
                        patch = obj.InterfacesRotor(iInt).patch1;
                        side = obj.InterfacesRotor(iInt).side1;
                        obj.intplotRotor{iInt} = plot(PlotPointsRotor{patch}{1, side+1}, PlotPointsRotor{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                    end
                    % Stator
                    for iMat = 1:numel(obj.Materials)
                        for patch = obj.Materials(iMat).PatchesStator
                            obj.srfplotStator{patch}=surf(PlotPointsStator{patch}{1, 1}, PlotPointsStator{patch}{2, 1}, zeros(size(obj.GeoPlotPointsStator{patch}{1, 1})), ...
                                "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.PlotFaceAlpha);
                        end
                    end
                    for iBnd = 1:numel(obj.BoundariesStator)
                        patch = obj.BoundariesStator(iBnd).patches;
                        side = obj.BoundariesStator(iBnd).faces;
                        obj.lplotStator{iBnd} = plot(PlotPointsStator{patch}{1, side+1}, PlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                    end
                    for iInt = 1:numel(obj.InterfacesStator)
                        patch = obj.InterfacesStator(iInt).patch1;
                        side = obj.InterfacesStator(iInt).side1;
                        obj.intplotStator{iInt} = plot(PlotPointsStator{patch}{1, side+1}, PlotPointsStator{patch}{2, side+1}, "Color", obj.PlotColor.Lines);
                    end

                case {'iter', 'done'}
                    % Rotor
                    for iMat = 1:numel(obj.Materials)
                        for patch = obj.Materials(iMat).PatchesRotor
                            set(obj.srfplotRotor{patch}, 'XData', PlotPointsRotor{patch}{1, 1}, 'YData', PlotPointsRotor{patch}{2, 1});
                        end
                    end
                    for iBnd = 1:numel(obj.BoundariesRotor)
                        patch = obj.BoundariesRotor(iBnd).patches;
                        side = obj.BoundariesRotor(iBnd).faces;
                        set(obj.lplotRotor{iBnd}, 'XData', PlotPointsRotor{patch}{1, side+1}, 'YData', PlotPointsRotor{patch}{2, side+1});
                    end
                    for iInt = 1:numel(obj.InterfacesRotor)
                        patch = obj.InterfacesRotor(iInt).patch1;
                        side = obj.InterfacesRotor(iInt).side1;
                        set(obj.intplotRotor{iInt}, 'XData', PlotPointsRotor{patch}{1, side+1}, 'YData', PlotPointsRotor{patch}{2, side+1});
                    end
                    % Stator
                    for iMat = 1:numel(obj.Materials)
                        for patch = obj.Materials(iMat).PatchesStator
                            set(obj.srfplotStator{patch}, 'XData', PlotPointsStator{patch}{1, 1}, 'YData', PlotPointsStator{patch}{2, 1});
                        end
                    end
                    for iBnd = 1:numel(obj.BoundariesStator)
                        patch = obj.BoundariesStator(iBnd).patches;
                        side = obj.BoundariesStator(iBnd).faces;
                        set(obj.lplotStator{iBnd}, 'XData', PlotPointsStator{patch}{1, side+1}, 'YData', PlotPointsStator{patch}{2, side+1});
                    end
                    for iInt = 1:numel(obj.InterfacesStator)
                        patch = obj.InterfacesStator(iInt).patch1;
                        side = obj.InterfacesStator(iInt).side1;
                        set(obj.intplotStator{iInt}, 'XData', PlotPointsStator{patch}{1, side+1}, 'YData', PlotPointsStator{patch}{2, side+1});
                    end
                otherwise
                    disp("Unknown state: "+string(state));
                    return
            end

            %warning('off', 'MATLAB:contour:ConstantData');


            xlim(obj.XLimits);
            ylim(obj.YLimits)
            title("Iteration: "+string(optimValues.iteration)+" Optimization Function: "+string(optimValues.fval));
            drawnow
            write2gif(obj.OptimizationFigure, optimValues.iteration+1, obj.SaveDirectory+"/Geometry_progress.gif");

        end

        function stop = storeOutput(obj, x, optimValues, state)
            x = x.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            for iParam = 1:obj.NumOptParameters
                if obj.OptimizationParameters(iParam).Name == "LENGTH"
                    length = x(iParam);
                end
                if obj.OptimizationParameters(iParam).Name == "ScaleR"
                    rscaling = x(iParam);
                end
            end
            disp(['Axial length: ' num2str(length) ' Radial scaling: ' num2str(rscaling)]);


            obj.OptimizationHistory(optimValues.funccount).iteration = optimValues.iteration;
            obj.OptimizationHistory(optimValues.funccount).fvalUsed = 1;

            stop = false;

            switch state
                case 'init'
                    
                case {'iter', 'done'}
                    opts = obj.Options;
                    opts.draw_geometry = false;
                    for iParam  = 1:numel (obj.OptimizationParameters)
                        opts.(obj.OptimizationParameters(iParam).Name) = x(iParam);
                        obj.OptimizationParameters(iParam).OptVal = x(iParam);
                    end
                    [srfRotor, ~, ~] = obj.RotorGeometryFunction(opts);
                    % Extract and update control points and surface
                    ControlPointsInit = obj.ControlPointsRotor;
                    for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                        ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                        ControlPointsInit(ind_loc, 1) = reshape(srfRotor(i).coefs(1, :, :, :)./srfRotor(i).coefs(4, :, :, :), [], 1);
                        ControlPointsInit(ind_loc, 2) = reshape(srfRotor(i).coefs(2, :, :, :)./srfRotor(i).coefs(4, :, :, :), [], 1);
                        ControlPointsInit(ind_loc, 3) = reshape(srfRotor(i).coefs(3, :, :, :)./srfRotor(i).coefs(4, :, :, :), [], 1);
                        ControlPointsInit(ind_loc, 4) = reshape(srfRotor(i).coefs(4, :, :, :), [], 1);
                    end
                    % Control point offsets
                    vals = x(obj.NumOptParameters+(1:obj.NumOptControlPointsRotor));
                    for iCtrPt = 1:numel(vals)
                        valsx = vals(iCtrPt).*cos(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                        valsy = vals(iCtrPt).*sin(obj.OptimizationControlPointsRotor(iCtrPt).Angle);
                        numbers = obj.OptimizationControlPointsRotor(iCtrPt).Number;
                        ControlPointsInit(numbers, 1) = ControlPointsInit(numbers, 1) + valsx;
                        ControlPointsInit(numbers, 2) = ControlPointsInit(numbers, 2) + valsy;
                    end
                    % update initial surface
                    for i  = 1:obj.MultipatchSpaceRotorGeo.npatch
                        ind_loc = obj.MultipatchSpaceRotorGeo.gnum{i};
                        srfRotor(i).coefs(1, :, :, :) = reshape(ControlPointsInit(ind_loc, 1).*ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                        srfRotor(i).coefs(2, :, :, :) = reshape(ControlPointsInit(ind_loc, 2).*ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                        srfRotor(i).coefs(3, :, :, :) = reshape(ControlPointsInit(ind_loc, 3).*ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                        srfRotor(i).coefs(4, :, :, :) = reshape(ControlPointsInit(ind_loc, 4), srfRotor(i).number);
                    end
                    nrbexport(srfRotor, obj.SaveDirectory + "/RotorOpt"+string(optimValues.iteration)+".txt")
                    [srfStator, ~, ~] = obj.StatorGeometryFunction(opts);
                    % Extract and update control points and surface
                    ControlPointsInit = obj.ControlPointsStator;
                    for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                        ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                        ControlPointsInit(ind_loc, 1) = reshape(srfStator(i).coefs(1, :, :, :)./srfStator(i).coefs(4, :, :, :), [], 1);
                        ControlPointsInit(ind_loc, 2) = reshape(srfStator(i).coefs(2, :, :, :)./srfStator(i).coefs(4, :, :, :), [], 1);
                        ControlPointsInit(ind_loc, 3) = reshape(srfStator(i).coefs(3, :, :, :)./srfStator(i).coefs(4, :, :, :), [], 1);
                        ControlPointsInit(ind_loc, 4) = reshape(srfStator(i).coefs(4, :, :, :), [], 1);
                    end
                    % Control point offsets
                    vals = x(obj.NumOptParameters+obj.NumOptControlPointsRotor+(1:obj.NumOptControlPointsStator));
                    for iCtrPt = 1:numel(vals)
                        valsx = vals(iCtrPt).*cos(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                        valsy = vals(iCtrPt).*sin(obj.OptimizationControlPointsStator(iCtrPt).Angle);
                        numbers = obj.OptimizationControlPointsStator(iCtrPt).Number;
                        ControlPointsInit(numbers, 1) = ControlPointsInit(numbers, 1) + valsx;
                        ControlPointsInit(numbers, 2) = ControlPointsInit(numbers, 2) + valsy;
                    end
                    % update initial surface
                    for i  = 1:obj.MultipatchSpaceStatorGeo.npatch
                        ind_loc = obj.MultipatchSpaceStatorGeo.gnum{i};
                        srfStator(i).coefs(1, :, :, :) = reshape(ControlPointsInit(ind_loc, 1).*ControlPointsInit(ind_loc, 4), srfStator(i).number);
                        srfStator(i).coefs(2, :, :, :) = reshape(ControlPointsInit(ind_loc, 2).*ControlPointsInit(ind_loc, 4), srfStator(i).number);
                        srfStator(i).coefs(3, :, :, :) = reshape(ControlPointsInit(ind_loc, 3).*ControlPointsInit(ind_loc, 4), srfStator(i).number);
                        srfStator(i).coefs(4, :, :, :) = reshape(ControlPointsInit(ind_loc, 4), srfStator(i).number);
                    end
                    nrbexport(srfStator, obj.SaveDirectory + "/StatorOpt"+string(optimValues.iteration)+".txt")
                    obj.OptimizationHistory(optimValues.funccount).Parameters = obj.OptimizationParameters;
                    OptimizationHistory = obj.OptimizationHistory;
                    save(obj.SaveDirectory+"/OptimizationHistory.mat", 'OptimizationHistory');            
                otherwise
                    disp("Unknown state: "+string(state));
            end
        end
        
        %% Optimization Constraints: radius of soft iron and monotonic incrasing B-Values
        function setConstraintFunction(obj, func)
            obj.ConstraintFunction = func;
        end

        function setRotorGeometryFunction(obj, func)
            obj.RotorGeometryFunction = func;
        end
        function setStatorGeometryFunction(obj, func)
            obj.StatorGeometryFunction = func;
        end
        
        function [c, ceq, gradc, gradceq] = optimizationConstraints(obj, x)
            [Tmean, dTmean] = obj.torqueConstraint(x);

            x = x.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            opts = obj.Options;
            opts.draw_geometry = false;
            for iParam  = 1:obj.NumOptParameters
                opts.(obj.OptimizationParameters(iParam).Name) = x(iParam);
            end
            c = obj.ConstraintFunction(opts);

            gradc = zeros(obj.NumOptParameters+obj.NumOptControlPointsRotor+obj.NumOptControlPointsStator, numel(c));
            
            for iParam = 1:obj.NumOptParameters
                opts1 = opts;
                Step = (obj.OptimizationParameters(iParam).MaxVal - obj.OptimizationParameters(iParam).MinVal)/1e3;
                opts1.(obj.OptimizationParameters(iParam).Name) = opts1.(obj.OptimizationParameters(iParam).Name) + Step;
                c1 = obj.ConstraintFunction(opts1);
                gradc(iParam, :) = reshape((c1-c)./Step, [], 1);
            end
            
            gradc = gradc.*(obj.UpperBounds-obj.LowerBounds);

            c = [c, Tmean];
            gradc = [gradc, dTmean];
            ceq = [];
            gradceq = [];
        end

        function plotOptimizationProgress(obj)
            Iter = [];
            MagnetMass = [];
            Tmean = [];
            Tsigma = [];
            for iIter = 1:numel(obj.OptimizationHistory)
                if obj.OptimizationHistory(iIter).iteration >= 0
                    Iter = [Iter, obj.OptimizationHistory(iIter).iteration];
                    MagnetMass = [MagnetMass, obj.OptimizationHistory(iIter).MagnetMass];
                    Tmean = [Tmean, obj.OptimizationHistory(iIter).Tmean];
                    Tsigma = [Tsigma, obj.OptimizationHistory(iIter).Tsigma];
                end
            end
            set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
            set(groot, 'defaultLegendInterpreter', 'latex');
            f = figure(42);
            clf
            set(gca, "FontSize", 24);
            Lw = 1.5;
            hold on
            grid on
%             yyaxis left
            plot(Iter, Tmean, 'LineStyle', '-', 'Color', TUDa_getColor("3d"), 'LineWidth', Lw)
            plot(Iter, MagnetMass*obj.MultCost, 'LineStyle', '-', 'Color', TUDa_getColor("5b"), 'LineWidth', Lw)
            plot(Iter, Tsigma*obj.MultTsig, 'LineStyle', '-', 'Color', TUDa_getColor("8a"), 'LineWidth', Lw)

            scatter(Iter, Tmean, 'filled', 'MarkerFaceColor', TUDa_getColor("3d"))
            scatter(Iter, MagnetMass*obj.MultCost, 'filled', 'MarkerFaceColor', TUDa_getColor("5b"))
            scatter(Iter, Tsigma*obj.MultTsig, 'filled', 'MarkerFaceColor', TUDa_getColor("8a"))
            legend(["Mean Torque", "Scaled Magnet Mass", "Scaled Torque std"])
            xlabel("Iteration", Interpreter="latex")
            ylim([0, 3]);
            xlim([0, max(Iter)])
            saveas(gca, obj.SaveDirectory +"/OptimizationProgress.fig")
            matlab2tikz('figurehandle', f, [char(obj.SaveDirectory) '/OptimizationProgress.tikz'])
        end

        function plotOptimizationTorque(obj, steps)
            set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
            set(groot, 'defaultLegendInterpreter', 'latex');
            figure(43)
%             clf
            set(gca, "FontSize", 20);
            hold on
            grid on
            colors = ["1a", "3a", "5a", "7a", "9a", "11a", "1c", "3c", "5c"];
            legendString = [];
            c = 1;
            Lw = 1.5;
            for iHist = 1:numel(obj.OptimizationHistory)
                if isempty(obj.OptimizationHistory(iHist).iteration)
                    continue
                end
                if any(obj.OptimizationHistory(iHist).iteration == steps)
                    plot(obj.OptimizationHistory(iHist).Angles, obj.OptimizationHistory(iHist).Torques, Color=TUDa_getColor(colors(c)), LineWidth=Lw);
                    legendString = [legendString, "Iteration "+ num2str(obj.OptimizationHistory(iHist).iteration)];
                    c = c +1;
                end
            end
            legend(legendString, 'Location', 'best');
            xlabel("Rotation angle", "Interpreter", "latex")
            ylabel("Torque", "Interpreter", "latex")
            %ylim([1.2, 1.85]);
            xlim([min(obj.OptimizationAngles), max(obj.OptimizationAngles)])
            saveas(gca, obj.SaveDirectory +"/OptimizationTorques.fig")
            matlab2tikz([char(obj.SaveDirectory) '/OptimizationTorques.tikz'])
        end

        function plotControlPointsUsed(obj)
            numbers = [];
            for iCtrlPts = 1:numel(obj.OptimizationControlPointsRotor)
                numbers = [numbers, reshape(obj.OptimizationControlPointsRotor(iCtrlPts).Number, 1, [])];
            end
            scatter(obj.ControlPointsRotor(numbers, 1), obj.ControlPointsRotor(numbers, 2), "filled", "red")

            numbers = [];
            for iCtrlPts = 1:numel(obj.OptimizationControlPointsStator)
                numbers = [numbers, reshape(obj.OptimizationControlPointsStator(iCtrlPts).Number, 1, [])];
            end
            scatter(obj.ControlPointsStator(numbers, 1), obj.ControlPointsStator(numbers, 2), "filled", "red")
%             text(obj.ControlPointsRotor(:, 1), obj.ControlPointsRotor(:, 2), string(1:size(obj.ControlPointsRotor, 1)));
        end

        function plotOptimizationHistoryTEMP(obj)
            %obj.FigureGeometry = figure('Name', 'Geometry Plot', 'NumberTitle', 'off', 'Position', [100 100 800 500]);
            obj.FigureGeometry = figure();
            obj.setRotationAngle(0)
            for iter = 0:200
                if ~isfile([obj.SaveDirectory + "/RotorOpt" + num2str(iter) + ".txt"])
                    break
                end
                clf
                title("Iteration "+string(iter))
                obj.GeometryRotor = mp_geo_load(char(obj.SaveDirectory + "/RotorOpt" + num2str(iter) + ".txt"));
                obj.CalculateGeoPlotPoints();
                hold on
                obj.plotPatches();
                obj.plotMaterialBoundaries();
                RI = find(vertcat(obj.OptimizationHistory.iteration) == iter);
                MA = obj.OptimizationHistory(RI).Parameters(4).OptVal;
                MA1 = pi/2-(deg2rad(MA/2) - pi/4);
                MA2 = MA1 -pi + deg2rad(MA);
                obj.addRemanenceRotor([9], MA1);
                obj.addRemanenceRotor([27], MA2);
                obj.plotRemanenceQuiver();
                axis equal
                axis off
                xlim(obj.XLimits);
                ylim(obj.YLimits);
                drawnow
%                 write2gif(obj.FigureGeometry, iter+1, obj.SaveDirectory+"/OptiHistory.gif");
                exportgraphics(obj.FigureGeometry, ['Progress' pad(num2str(iter), 3, 'left', '0') '.png']);
            end
        end

        function loadOptimizationHisotry(obj)
            [file, path] = uigetfile();
            obj.SaveDirectory = path;
            obj.OptimizationHistory = load([path file]).OptimizationHistory;
        end

        %% srf calculation for optimization
        function GeoPlotPoints = CalculateGeoPlotPointsSrf(obj, srf)
            for iPatch = 1:numel(srf)
%                 ParamPlotPoints{i} = {linspace(0, 1, srf(i).order(1)* srf(i).number(1)), ...
%                     linspace(0, 1, srf(i).order(2)* srf(i).number(2))};
                pointsX = unique(kntrefine(srf(iPatch).knots{1}, obj.PlotResolutionX, srf(iPatch).order(1), srf(iPatch).order(1)-1));
                pointsY = unique(kntrefine(srf(iPatch).knots{2}, obj.PlotResolutionY, srf(iPatch).order(2), srf(iPatch).order(2)-1));

                ParamPlotPoints{iPatch} = {pointsX , pointsY};

                GeoPlotPoints{iPatch} = obj.getGeoPlotPointsSrf(srf(iPatch), ParamPlotPoints{iPatch}, 0);
            end
        end

        function pnts = getGeoPlotPointsSrf(obj, nurbs, pts, rot)
            xpts = pts{1};%linspace(0, 1, geo.nurbs.order(1)* geo.nurbs.number(1));
            ypts = pts{2};%linspace(0, 1, geo.nurbs.order(2)* geo.nurbs.number(2));
        
            F_geo = nrbeval(nurbs, {xpts, ypts});
            X = squeeze(F_geo(1, :, :))*cos(rot) - squeeze(F_geo(2, :, :))*sin(rot);
            Y = squeeze(F_geo(2, :, :))*cos(rot) + squeeze(F_geo(1, :, :))*sin(rot);
            pnts{1, 1} = X; 
            pnts{2, 1} = Y; 
        
            bnd = nrbeval(nrbextract(nurbs, 1), ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 2} = X;
            pnts{2, 2} = Y;
        
            bnd = nrbeval(nrbextract(nurbs, 2), ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 3} = X;
            pnts{2, 3} = Y;
        
            bnd = nrbeval(nrbextract(nurbs, 3), xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 4} = X;
            pnts{2, 4} = Y;
        
            bnd = nrbeval(nrbextract(nurbs, 4), xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 5} = X;
            pnts{2, 5} = Y;  
        end
    end
end
