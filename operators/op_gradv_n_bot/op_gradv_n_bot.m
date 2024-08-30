% OP_GRADV_N_BOT: assemble the right-hand side vector r = [r(i)], with  r(i) = (grad(v)_i, n_bot(alpha)).
% where n_bot is the perpendicular normal vector to the magnetization
% direction given by alpha: n_bot = [-sin(alpha); cos(alpha)]
%   rhs = op_gradv_n_bot (spv, msh, alpha);
% INPUT:

%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   alpha: magnetization angle in rad
%
% OUTPUT:
%
%   rhs:    assembled right hand side

function rhs = op_gradv_n_bot (spv, msh, alpha)

    gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
    rhs   = zeros (spv.ndof, 1);
    
    for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
     jacdet_weights = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel), 1, msh.nqn);
    
     gradv_iel = reshape (gradv(:, :, :, :, iel), spv.ncomp*msh.ndim, msh.nqn, spv.nsh_max, 1);
    
     aux_val = gradv_iel .* [-sin(alpha); cos(alpha)] .*jacdet_weights;
    
     rhs_loc = sum (sum (aux_val, 1), 2);
    
     indices = find (spv.connectivity(:,iel));
     rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
     rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:); 
    else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradv_n_bot: singular map in element number %d', iel)
    end
    end

end