function T = op_Br_Bt_dGamma_tp(space, msh, u)

    for idim = 1:msh.ndim
        size1 = size (space.sp_univ(idim).connectivity);
        if (size1(2) ~= msh.nel_dir(idim))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    T = 0;
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp_col = sp_evaluate_col(space, msh_col, 'value', true, 'gradient', true);
        T = T + op_Br_Bt_dGamma(sp_col, msh_col, u);
    end
end