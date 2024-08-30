% OP_D_DC_OMEGA: assemble the derivative of the patch volume w.r.t the control points of g.
%
%   mat = op_D_DC_Omega (spg, msh);
%   [row, col, values] = op_D_DC_Omega (spg, msh);
%
% INPUT:
%
%   spg:   structure representing the PARAMETRIC space of geometry functions (see sp_scalar/sp_evaluate_col_param/sp_precompute_param)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   row:    nonzero row entries
%   col:    nonzero column entries
%   values: values of the nonzero entries

function varargout = op_D_DC_Omega(spgP, msh)
    
    gradG = reshape (spgP.shape_function_gradients, spgP.ncomp, [], msh.nqn, spgP.nsh_max, msh.nel); % Parametric
    
    rows = zeros (msh.nel * spgP.nsh_max * msh.ndim, 1);
    dims = zeros (msh.nel * spgP.nsh_max * msh.ndim, 1);
    values = zeros (msh.nel * spgP.nsh_max * msh.ndim, 1);

    jacdet_weights = msh.jacdet .* msh.quad_weights;
    
    ncounter = 0;
    for iel = 1:msh.nel
        if (all (msh.jacdet(:, iel)))

            gradG_iel_mod = reshape (gradG(:, :, :, :, iel), msh.ndim, 1, msh.nqn, spgP.nsh_max);
            jacdet_iel = reshape (jacdet_weights(:, iel), [1, msh.nqn, 1, 1]);
            
            Jinv = pageinv(msh.geo_map_jac(:, :, :, iel)); 

            TraceJinvTG = reshape(sum(Jinv .* gradG_iel_mod, 1), msh.ndim, msh.nqn, spgP.nsh_max);
            dVdP = jacdet_iel .* TraceJinvTG;

            elementary_values = reshape(sum(dVdP, 2), [msh.ndim, spgP.nsh_max]);
        
            [dims_loc, rows_loc] = ndgrid (1:msh.ndim, spgP.connectivity(:, iel));
            indices = dims_loc & rows_loc;

            rows(ncounter+(1:spgP.nsh(iel)*msh.ndim)) = rows_loc(indices);
            dims(ncounter+(1:spgP.nsh(iel)*msh.ndim)) = dims_loc(indices);
            
            values(ncounter+(1:spgP.nsh(iel)*msh.ndim)) = elementary_values(indices);
            ncounter = ncounter + spgP.nsh(iel)*msh.ndim;
    
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_D_DC_Omega: singular map in element number %d', iel)
        end
    end
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows(1:ncounter), dims(1:ncounter),...
            values(1:ncounter), spgP.ndof, msh.ndim);
    elseif (nargout == 3)
        varargout{1} = rows(1:ncounter);
        varargout{2} = dims(1:ncounter);
        varargout{3} = values(1:ncounter);
    else
        error ('op_D_DC_gradu_gradv: wrong number of output arguments')
    end

end

