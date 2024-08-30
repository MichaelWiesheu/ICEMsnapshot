% Faster version of op_f_v with multiple functions for f

function varargout = op_fs_v (spv, msh, fs)

    for idim = 1:msh.rdim
        x{idim} = reshape (msh.geo_map(idim, :, :), msh.nqn, msh.nel);
    end    
    fev = fs (x{:});
    nfs = size(fev, 1);
    
    shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
    
    rows = zeros (msh.nel * nfs * spv.nsh_max, 1);
    cols = zeros (msh.nel * nfs * spv.nsh_max, 1);
    values = zeros (msh.nel * nfs * spv.nsh_max, 1);
    
    jacdet_weights = msh.jacdet .* msh.quad_weights;
    
    ncounter = 0;
    for iel = 1:msh.nel
        if (all (msh.jacdet(:, iel)))
            fevel = reshape (fev(:, :, iel), nfs, msh.nqn);
            shpv_iel = reshape (shpv(:, :, :, iel), spv.ncomp, msh.nqn, spv.nsh_max);
    
            jacdet_iel = reshape (jacdet_weights(:, iel), [1, msh.nqn, 1]);
    
            elementary_values = reshape (sum (shpv_iel.*fevel.*jacdet_iel, 2), nfs, spv.nsh_max);
    
            [cols_loc, rows_loc] = ndgrid (1:nfs, spv.connectivity(:, iel));
            indices =  cols_loc & rows_loc;
            rows(ncounter+(1:nfs*spv.nsh(iel))) = rows_loc(indices);
            cols(ncounter+(1:nfs*spv.nsh(iel))) = cols_loc(indices);
            values(ncounter+(1:nfs*spv.nsh(iel))) = elementary_values(indices);
            ncounter = ncounter + nfs*spv.nsh(iel);
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_fs_v: singular map in element number %d', iel)
        end
    end
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
            values(1:ncounter), spv.ndof, nfs);
    elseif (nargout == 3)
        varargout{1} = rows(1:ncounter);
        varargout{2} = cols(1:ncounter);
        varargout{3} = values(1:ncounter);
    else
        error ('op_fs_v: wrong number of output arguments')
    end

end


