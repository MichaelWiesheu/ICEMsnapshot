% OP_GRADV_N_BOT_MP: assemble the right-hand side vector r = [r(i)], with  r(i) = (grad(v)_i, n_bot(alpha)).
% where n_bot is the perpendicular normal vector to the magnetization
% direction given by alpha: n_bot = [-sin(alpha); cos(alpha)] for a
% multipatch domain
%   rhs = op_gradv_n_bot_mp (spv, msh, alpha);
% INPUT:

%   spv:     object representing the function space (see sp_multipatch)
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   alpha:   magnetization angle in rad
%
% OUTPUT:
%
%   rhs:    assembled right hand side

function rhs = op_gradv_n_bot_mp (space, msh, alpha, patch_list)

  if (nargin < 4)
    patch_list = 1:msh.npatch;
  end

  if (space.npatch ~= msh.npatch)
    error ('op_gradv_n_bot_mp: the number of patches does not coincide')
  end
  
  rhs = zeros (space.ndof, 1);
  for iptc = patch_list
    rhs_loc = op_gradv_n_bot_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, alpha);
    
    if (~isempty (space.dofs_ornt))
      rhs_loc = space.dofs_ornt{iptc}(:) .* rhs_loc(:);
    end
    rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
  end

end