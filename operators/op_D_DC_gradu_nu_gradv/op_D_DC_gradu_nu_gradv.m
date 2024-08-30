% OP_D_DC_GRADU_GRADV: assemble the derivative of the stiffness matrix A = [a(i,j)], a(i,j) = (nu grad u_j, grad v_i) w.r.t the control points of g.
%
%   mat = op_gradu_gradv (spu, spv, spg, msh, u, mat);
%   [indices, values] = op_gradu_gradv (spu, spv, spg, msh, u, mat);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   spg:   structure representing the space of geometry functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   u:     discrete magnetic vector potential
%   mat:   material class with reluctivity functions
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   indices:  indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_D_DC_gradu_nu_gradv(spu, spv, spgP, msh, u, mat)
    
    gradG = reshape (spgP.shape_function_gradients, spgP.ncomp, [], msh.nqn, spgP.nsh_max, msh.nel); % Parametric
    
    rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    tens = zeros (msh.nel * spu.nsh_max * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    dims = zeros (msh.nel * spu.nsh_max * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    values = zeros (msh.nel * spu.nsh_max * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    
    ncounter = 0;
    for iel = 1:msh.nel
        if (all (msh.jacdet(:, iel)))
            gradu_iel = reshape (spu.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, 1, spu.nsh_max);
            gradv_iel = reshape (spv.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, spv.nsh_max, 1);
            %gradG_iel = reshape (gradG(:,:,:,:,iel), msh.ndim, msh.nqn, 1, 1, spgP.nsh_max);

            gradG_iel_mod = reshape (gradG(:,:,:,:,iel), msh.ndim, 1, msh.nqn, 1, 1, spgP.nsh_max);
            jacdet_iel = reshape (msh.jacdet(:,iel), [1, msh.nqn, 1, 1]);
            weight_iel = reshape (msh.quad_weights(:,iel), [1, msh.nqn, 1, 1]);

            u_iel = reshape(u(spu.connectivity(:,iel)), [1, 1, 1, spu.nsh_max]);
            gradNi_ui = sum(gradu_iel.*u_iel, 4);
            
            B_iel = (sum(gradNi_ui.^2, 1)).^0.5;
            nu_iel = mat.getNuNonlinear(B_iel);
            dnudB_iel = mat.getNuPrimeNonlinear(B_iel);
            
            Jinv = pageinv(msh.geo_map_jac(:, :, :, iel)); 
            
            JinvTgradG = reshape(sum(Jinv.*gradG_iel_mod, 1), [msh.ndim, msh.nqn, 1, 1, spgP.nsh_max]);

            dK1dP = -weight_iel.*jacdet_iel.*nu_iel.*sum(JinvTgradG.*gradv_iel, 1).*gradu_iel;
            dK2dP = -weight_iel.*jacdet_iel.*nu_iel.*sum(JinvTgradG.*gradu_iel, 1).*gradv_iel;

            TraceJinvTG = reshape(sum(Jinv.*gradG_iel_mod, 1), msh.ndim, msh.nqn, 1, 1, spgP.nsh_max);
            dK3dP = weight_iel.*jacdet_iel.*nu_iel.*TraceJinvTG.*sum(gradu_iel.*gradv_iel, 1);

            dK4dP(1,:,:,:,:) = -weight_iel.*jacdet_iel.*sum(gradu_iel.*gradv_iel, 1).*dnudB_iel./B_iel.*sum(sum(gradu_iel.*u_iel,4).*(JinvTgradG.*sum(gradu_iel(1,:,:,:).*u_iel,4)),1);  
            dK4dP(2,:,:,:,:) = -weight_iel.*jacdet_iel.*sum(gradu_iel.*gradv_iel, 1).*dnudB_iel./B_iel.*sum(sum(gradu_iel.*u_iel,4).*(JinvTgradG.*sum(gradu_iel(2,:,:,:).*u_iel,4)),1); 

            dKdP = dK1dP + dK2dP + dK3dP + dK4dP;
            
            elementary_values = reshape(sum(dKdP, 2), [msh.ndim, spv.nsh_max, spu.nsh_max, spgP.nsh_max]);
        
            [dims_loc, rows_loc, cols_loc, tens_loc] = ndgrid (1:msh.ndim, spv.connectivity(:,iel), spu.connectivity(:,iel), spgP.connectivity(:,iel));
            indices = dims_loc & rows_loc & cols_loc & tens_loc;
            rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = rows_loc(indices);
            cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = cols_loc(indices);
            tens(ncounter+(1:spu.nsh(iel)*spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = tens_loc(indices);
            dims(ncounter+(1:spu.nsh(iel)*spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = dims_loc(indices);
            
            values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = elementary_values(indices);
            ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel)*spgP.nsh(iel)*msh.ndim;
    
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_D_DC_gradu_nu_gradv: singular map in element number %d', iel)
        end
    end
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows(1:ncounter), cols(1:ncounter), tens(1:ncounter), dims(1:ncounter)],...
            values(1:ncounter), [spu.ndof, spv.ndof, spgP.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows(1:ncounter), cols(1:ncounter), tens(1:ncounter), dims(1:ncounter)];
        varargout{2} = values(1:ncounter);
    else
        error ('op_D_DC_gradu_nu_gradv: wrong number of output arguments')
    end

end

