% Function for permanent magnet
% Int gradv [-sin(alpha); cos(alpha)] dOmega
% RHS of magnet system, where remanence direction is given by alpha

function rhs = op_d_dalpha_gradv_n_bot (spv, msh, alpha)

    gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
    rhs   = zeros (spv.ndof, 1);
    
    for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
     jacdet_weights = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel), 1, msh.nqn);
    
     gradv_iel = reshape (gradv(:,:,:,:,iel), spv.ncomp*msh.ndim, msh.nqn, spv.nsh_max, 1);
    
     aux_val = gradv_iel .* [-cos(alpha); -sin(alpha)] .*jacdet_weights;
    
     rhs_loc = sum (sum (aux_val, 1), 2);
    
     indices = find (spv.connectivity(:,iel));
     rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
     rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:); 
    else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_d_dalpha_gradv_n_bot: singular map in element number %d', iel)
    end
    end

end