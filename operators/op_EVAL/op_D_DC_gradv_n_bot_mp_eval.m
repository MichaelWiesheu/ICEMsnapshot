function A = op_D_DC_gradv_n_bot_mp_eval(spv, spvEval, spg, spgEval, msh, mshEval, alpha, patch_list)

  inds = cell(msh.npatch, 1);
  vals = cell(msh.npatch, 1);
  
  for iPatch = patch_list
    [ind, v] = op_D_DC_gradv_n_bot(spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, alpha);

    ind(:,1) = spv.gnum{iPatch}(ind(:,1));
    ind(:,2) = spg.gnum{iPatch}(ind(:,2));

%     if (~isempty (spv.dofs_ornt))
%       v = spv.dofs_ornt{iPatch}(ind(:,1))' .* v;
%     end

    inds{iPatch} = ind;
    vals{iPatch} = v;

  end

  A = sptensor(cell2mat(inds),cell2mat(vals), [spv.ndof,spg.ndof,msh.ndim]);

end