% TBD!!
% % % OP_D_DC_GRADU_GRADV_MP: assemble the derivative of the stiffness matrix A_dP = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), in a multipatch domain with respect to the geometry control points of g_k.
% % %
% % %   mat = op_gradu_gradv_mp (spu, spv, spg, msh, [epsilon], [patches]);
% % %
% % % INPUT:
% % %
% % %   spu:     object representing the space of trial functions (see sp_multipatch)
% % %   spv:     object representing the space of test functions (see sp_multipatch)
% % %   spg:     object representing the geometry space
% % %   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
% % %   epsilon: function handle to compute the diffusion coefficient. Equal to one if left empty.
% % %   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
% % %
% % % OUTPUT:
% % %
% % %   mat:    derivative of assembled stiffness matrix by control points

function A = op_D_DC_gradv_n_bot_mp(spv, spg, msh, alpha, patch_list)

    if (nargin < 5)
        patch_list = 1:msh.npatch;
    end
    
    if ((spv.npatch~= spg.npatch)  || (spg.npatch ~= msh.npatch))
        error ('op_D_DC_gradv_n_bot_mp: the number of patches does not coincide')
    end
    
    rws = [];
    cls = [];
    dms = [];
    vals = [];
    ncounter = 0;
    for iptc = patch_list
    
        [indices, vs] = op_D_DC_gradv_n_bot_tp (spv.sp_patch{iptc}, spg.sp_patch{iptc}, msh.msh_patch{iptc}, alpha);

        if isempty(indices)
            continue
        end
    
        rs = indices(:,1);
        cs = indices(:,2);
        dm = indices(:,3);
    
        rws(ncounter+(1:numel (rs))) = spv.gnum{iptc}(rs);
        cls(ncounter+(1:numel (cs))) = spg.gnum{iptc}(cs);
        dms(ncounter+(1:numel (dm))) = dm;
    
        if (~isempty (spv.dofs_ornt))
            vs = vs .* spv.dofs_ornt{iptc}(rs)';
            warning("op_D_DC_gradv_n_bot_mp: dofs_ornt not tested yet")
        end
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_gradv_n_bot_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iptc}(cs)';
        end
    
        vals(ncounter+(1:numel (rs))) = vs;
        ncounter = ncounter + numel (rs);
    end

    rws = reshape(rws, [], 1);
    cls = reshape(cls, [], 1);
    dms = reshape(dms, [], 1);
    vals = reshape(vals, [], 1);
    indices = [rws, cls, dms];
    
    A = sptensor (indices, vals, [spv.ndof, spg.ndof, msh.ndim]);

end
