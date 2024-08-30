
function rhs = op_gradv_n_bot_mp_eval(spv, spvEval, mshEval, alpha, patch_list)

  rhs = zeros (spv.ndof, 1);
  
  for iPatch = patch_list
    rhs_loc = op_gradv_n_bot(spvEval{iPatch}, mshEval{iPatch}, alpha);

    if (~isempty (spv.dofs_ornt))
      rhs_loc = spv.dofs_ornt{iptc}(:) .* rhs_loc(:);
    end

    rhs(spv.gnum{iPatch}) = rhs(spv.gnum{iPatch}) + rhs_loc;

  end
end