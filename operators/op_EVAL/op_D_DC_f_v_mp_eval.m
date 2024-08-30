function A = op_D_DC_f_v_mp_eval(spv, spvEval, spg, spgEval, msh, mshEval, coeff, patch_list)

  inds = cell(msh.npatch, 1);
  vals = cell(msh.npatch, 1);
  
  for iPatch = patch_list
    [ind, v] = op_D_DC_f_v(spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, coeff);

    ind(:,1) = spv.gnum{iPatch}(ind(:,1));
    ind(:,2) = spg.gnum{iPatch}(ind(:,2));

    if (~isempty (spv.dofs_ornt))
        warning("op_D_DC_f_v_mp_eval: dofs_ornt not tested!")
        vs = vs .* spv.dofs_ornt{iptc}(ind(:,1))';
    end

    inds{iPatch} = ind;
    vals{iPatch} = v;

  end

  A = sptensor(cell2mat(inds),cell2mat(vals), [spv.ndof,spg.ndof,msh.ndim]);
end