% TBD
% % % op_D_DC_gradv_n_bot_tp: assemble the derivative of the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), w.r.t the control points g_p exploiting the tensor product structure.
% % %
% % %   mat = op_gradu_gradv_tp (spu, spv, spg, msh, [epsilon]);
% % %   [indices, values] = op_gradu_gradv_tp (spu, spv, spg, msh, [epsilon]);
% % %
% % % INPUT:
% % %
% % %   spu:     object representing the space of trial functions (see sp_vector)
% % %   spv:     object representing the space of test functions (see sp_vector)
% % %   spg:     object representing the space of the geometry functions (see sp_scalar)
% % %   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
% % %   epsilon: function handle to compute the diffusion coefficient (optional)
% % %
% % % OUTPUT:
% % %
% % %   mat:      assembled stiffness tensor
% % %   indices:  indices of the nonzero entries
% % %   values:   values of the nonzero entries

function varargout = op_D_DC_gradv_n_bot_tp(space1, space2, msh, alpha)

  for idim = 1:msh.ndim
    size1 = size (space1.sp_univ(idim).connectivity);
    size2 = size (space2.sp_univ(idim).connectivity);
    if ((size1(2) ~= size2(2)) || (size2(2) ~= msh.nel_dir(idim)))
      error ('One of the discrete spaces is not associated to the mesh')
    end
  end

  A = sptensor([],[],[space1.ndof, space2.ndof, msh.ndim]);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col(space1, msh_col, 'value', true, 'gradient', true);
    sp2_col_param = sp_evaluate_col_param(space2, msh_col, 'value', true, 'gradient', true);

    A = A + op_D_DC_gradv_n_bot(sp1_col, sp2_col_param, msh_col, alpha);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 2)
    [indices, vals] = find (A);
    varargout{1} = indices;
    varargout{2} = vals;
  end

end