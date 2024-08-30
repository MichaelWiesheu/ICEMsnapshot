% OP_D_DC_GRADU_GRADV_MP: assemble the derivative of the stiffness matrix A_dP = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i), in a multipatch domain with respect to the geometry control points of g_k.
%
%   mat = op_gradu_gradv_mp (spu, spv, spg, msh, [epsilon], [patches]);
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_multipatch)
%   spv:     object representing the space of test functions (see sp_multipatch)
%   spg:     object representing the geometry space
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   epsilon: function handle to compute the diffusion coefficient. Equal to one if left empty.
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    derivative of assembled stiffness matrix by control points

function A = op_D_DC_gradu_gradv_mp(spu, spv, spg, msh, coeff, patch_list)

    if (nargin < 6)
        patch_list = 1:msh.npatch;
    end
    
    if ((spu.npatch~= spv.npatch)  || (spv.npatch~= spg.npatch)  || (spu.npatch ~= msh.npatch))
        error ('op_D_DC_gradu_gradv_mp: the number of patches does not coincide')
    end
    
    rws = [];
    cls = [];
    tns = [];
    dms = [];
    vals = [];
    ncounter = 0;
    for iptc = patch_list
    
        if (nargin < 5 || isempty (coeff))
          [indices, vs] = op_D_DC_gradu_gradv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, spg.sp_patch{iptc}, msh.msh_patch{iptc});
        else
          [indices, vs] = op_D_DC_gradu_gradv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, spg.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
        end

        if isempty(indices)
            continue
        end
    
        rs = indices(:,1);
        cs = indices(:,2);
        ts = indices(:,3);
        dm = indices(:,4);
    
        rws(ncounter+(1:numel (rs))) = spu.gnum{iptc}(rs);
        cls(ncounter+(1:numel (cs))) = spv.gnum{iptc}(cs);
        tns(ncounter+(1:numel (ts))) = spg.gnum{iptc}(ts);
        dms(ncounter+(1:numel (dm))) = dm;
    
        if (~isempty (spu.dofs_ornt))
            vs = spu.dofs_ornt{iptc}(rs)' .* vs;
        end
        if (~isempty (spv.dofs_ornt))
            vs = vs .* spv.dofs_ornt{iptc}(cs)';
        end
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_gradu_gradv_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iptc}(ts)';
        end
    
        vals(ncounter+(1:numel (rs))) = vs;
        ncounter = ncounter + numel (rs);
    end

    rws = reshape(rws, [], 1);
    cls = reshape(cls, [], 1);
    tns = reshape(tns, [], 1);
    dms = reshape(dms, [], 1);
    vals = reshape(vals, [], 1);
    indices = [rws, cls, tns, dms];
    
    A = sptensor (indices, vals, [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);

end
