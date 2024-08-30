

function varargout = op_D_DC_gradv_n_bot (spv, spgP, msh, alpha)

    gradG = reshape (spgP.shape_function_gradients, spgP.ncomp, [], msh.nqn, spgP.nsh_max, msh.nel); % Parametric

    rows = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    cols = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    dims = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    values = zeros (msh.nel * spv.nsh_max * spgP.nsh_max * msh.ndim, 1);
    
    ncounter = 0;
    for iel = 1:msh.nel
        if (all (msh.jacdet(:, iel)))
        
        gradv_iel = reshape (spv.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, spv.nsh_max, 1);
        gradG_iel = reshape (spgP.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, 1, spgP.nsh_max);
        %gradG_iel = reshape (gradG(:,:,:,:,iel), msh.ndim, msh.nqn, 1, 1, spgP.nsh_max);
        
        gradG_iel_mod = reshape (gradG(:,:,:,:,iel), msh.ndim, 1, msh.nqn, 1, spgP.nsh_max);
        
        jacdet_iel = reshape (msh.jacdet(:,iel), [1, msh.nqn, 1, 1]);
        weight_iel = reshape (msh.quad_weights(:,iel), [1, msh.nqn, 1, 1]);
        
        Jinv = pageinv(msh.geo_map_jac(:, :, :, iel)); 
        
        JinvTgradG = reshape(sum(Jinv.*gradG_iel_mod, 1), [msh.ndim, msh.nqn, 1, spgP.nsh_max]);
        df1dP = -weight_iel.*jacdet_iel.*sum(JinvTgradG.*[-sin(alpha); cos(alpha)], 1).*gradv_iel;

        TraceJinvTG = reshape(sum(Jinv.*gradG_iel_mod, 1), msh.ndim, msh.nqn, 1, spgP.nsh_max);
        df2dP = weight_iel.*jacdet_iel.*TraceJinvTG.*sum([-sin(alpha); cos(alpha)].*gradv_iel, 1);

%         test1 = -weight_iel.*jacdet_iel.*sum(JinvTgradG.*gradv_iel, 1).*[-sin(alpha); cos(alpha)];
%         test2 = weight_iel.*jacdet_iel.*TraceJinvTG.*sum([-sin(alpha); cos(alpha)].*gradv_iel, 1);
%         test3 = -weight_iel.*jacdet_iel.*sum(JinvTgradG.*[-sin(alpha); cos(alpha)], 1).*gradv_iel;

        dfdP = df1dP + df2dP;
%         dfdP = test2+test3;

        elementary_values = reshape(sum(dfdP, 2), [msh.ndim, spv.nsh_max, spgP.nsh_max]);
    
        [dims_loc, rows_loc, cols_loc] = ndgrid (1:msh.ndim, spv.connectivity(:,iel), spgP.connectivity(:,iel));
        indices = dims_loc & rows_loc & cols_loc;
        rows(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = rows_loc(indices);
        cols(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = cols_loc(indices);
        dims(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = dims_loc(indices);
        
        values(ncounter+(1:spv.nsh(iel)*spgP.nsh(iel)*msh.ndim)) = elementary_values(indices);
        ncounter = ncounter + spv.nsh(iel)*spgP.nsh(iel)*msh.ndim;
    
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_D_DC_gradv_n_bot: singular map in element number %d', iel)
        end
    end
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows(1:ncounter), cols(1:ncounter), dims(1:ncounter)],...
            values(1:ncounter), [spv.ndof, spgP.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows(1:ncounter), cols(1:ncounter), dims(1:ncounter)];
        varargout{2} = values(1:ncounter);
    else
        error ('op_D_DC_u_v: wrong number of output arguments')
    end

end