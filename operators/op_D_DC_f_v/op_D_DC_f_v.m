% OP_D_DC_F_V: assemble the derivative of the rhs b = (v_i) w.r.t the control points of g.
%
%
% INPUT:
%
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   spg:   structure representing the PARAMETRIC space of geometry functions (see sp_scalar/sp_evaluate_col_param/sp_precompute_param)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   epsilon: CONSTANT coefficient
%
% OUTPUT:
%
%   mat:    assembled derivative tensor
%   ind:    sptensor indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_D_DC_f_v(spv, spgP, msh, coeff)
    
    gradG = reshape (spgP.shape_function_gradients, spgP.ncomp, [], msh.nqn, spgP.nsh_max, msh.nel); % Parametric
    
    rows = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    cols = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    dims = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    values = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);

    jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
    
    ncounter = 0;
    for iel = 1:msh.nel
        if (all (msh.jacdet(:, iel)))
            shpv_iel = reshape (spv.shape_functions(:, :, iel), 1, msh.nqn, spv.nsh_max);
            gradv_iel = reshape (spv.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, spv.nsh_max, 1);
            gradG_iel = reshape (spgP.shape_function_gradients(:, :, :, iel), msh.ndim, 1, msh.nqn, spgP.nsh_max);

            G_iel = reshape(spgP.shape_functions(:,:,iel), 1, msh.nqn, 1, spgP.nsh_max);
            jacdet_iel = reshape (jacdet_weights(:, iel), [1, msh.nqn, 1, 1]);
            
            db1dC = jacdet_iel .* gradv_iel.*G_iel;
            
            Jinv = pageinv(msh.geo_map_jac(:, :, :, iel)); 
            TraceJinvTG = reshape(sum(Jinv .* gradG_iel, 1), msh.ndim, msh.nqn, 1, spgP.nsh_max);
            db2dC = jacdet_iel .* shpv_iel .* TraceJinvTG ;
            
            dbdC = db2dC; %db1dC + 
%             b1 = squeeze(sum(db1dC,2));
%             b2 = squeeze(sum(db2dC,2));


            elementary_values = reshape(sum(dbdC, 2), [msh.ndim, spv.nsh_max, spgP.nsh_max]);
        
            [dims_loc, rows_loc, cols_loc] = ndgrid (1:msh.ndim, spv.connectivity(:, iel), spgP.connectivity(:, iel));
            indices = dims_loc & rows_loc & cols_loc;
            rows(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = rows_loc(indices);
            cols(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = cols_loc(indices);
            dims(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = dims_loc(indices);
            
            values(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = elementary_values(indices);
            ncounter = ncounter + spv.nsh(iel)*spgP.nsh(iel)*msh.ndim;
    
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_D_DC_f_v: singular map in element number %d', iel)
        end
    end
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows(1:ncounter), cols(1:ncounter), dims(1:ncounter)],...
            values(1:ncounter), [spv.ndof, spgP.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows(1:ncounter), cols(1:ncounter), dims(1:ncounter)];
        varargout{2} = values(1:ncounter);
    else
        error ('op_D_DC_f_v: wrong number of output arguments')
    end

end

