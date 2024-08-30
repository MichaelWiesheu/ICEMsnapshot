% sp_eval1: spinoff sp_eval returning the shape functions at the mesh qp
% points
% tested only for value and gradients

function [shp] = shp_eval1 (space, geometry, npts, options)


  if (nargin < 4)
    options = {'value'};
  end
  if (~iscell (options))
    options = {options};
  end
  nopts = numel (options);

  ndim = numel (space.knots);
  
  endpoints = zeros (2, ndim);
  if (isfield (geometry, 'nurbs'))
    nurbs = geometry.nurbs;
    if (ndim == 1)
      nurbs.knots = {nurbs.knots};
    end
    for idim=1:ndim
      endpoints(:,idim) = nurbs.knots{idim}([nurbs.order(idim), end-nurbs.order(idim)+1]);
    end
    clear nurbs
  elseif (isfield (struct(space), 'knots'))
    degree = space.degree;
    for idim=1:ndim
      endpoints(:,idim) = space.knots{idim}([degree(idim)+1, end-degree(idim)]);
    end
  else
    endpoints(2,:) = 1;
  end
  
  if (iscell (npts))
    pts = npts;
    npts = cellfun (@numel, pts);
  elseif (isvector (npts))
    if (numel (npts) == 1)
      npts = npts * ones (1,ndim);
    end
    for idim = 1:ndim
      pts{idim} = linspace (endpoints(1,idim), endpoints(2,idim), npts(idim));
    end
  end

  for jj = 1:ndim
    pts{jj} = pts{jj}(:)';
    if (numel (pts{jj}) > 1)
      brk{jj} = [endpoints(1,jj), pts{jj}(1:end-1) + diff(pts{jj})/2, endpoints(2,jj)];
    else
      brk{jj} = endpoints(:,jj).';
    end
  end

  msh = msh_cartesian (brk, pts, [], geometry, 'boundary', false);
  sp  = space.constructor (msh);

  
  value = false; grad = false; laplacian = false; hessian = false;
  
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        shp{iopt} = zeros (numel(msh.qn{1}), numel(msh.qn{2}), sp.ndof);
        value = true;

      case 'gradient'
        shp{iopt} = zeros (2, numel(msh.qn{1}), numel(msh.qn{2}), sp.ndof);
        grad = true;
        
      case 'laplacian'
        shp{iopt} = zeros (msh.nqn, sp.ndof);
        laplacian = true;

      case 'hessian'
        shp{iopt} = zeros (msh.rdim, msh.rdim, msh.nqn, sp.ndof);
        hessian = true;
    end
  end
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (sp, msh_col, 'value', value, 'gradient', grad, ...
          'laplacian', laplacian, 'hessian', hessian);  
    for iopt = 1:nopts
        switch (lower (options{iopt}))
            case 'value'
                for jel = 1:msh.nel_dir(2)
                    shp{iopt}(iel, jel, sp_col.connectivity(:,jel)) = sp_col.shape_functions(1,:,jel); 
                end
            case 'gradient'
                for jel = 1:msh.nel_dir(2)
                    shp{iopt}(1, iel, jel,sp_col.connectivity(:,jel)) = sp_col.shape_function_gradients(1,1,:,jel);
                    shp{iopt}(2, iel, jel,sp_col.connectivity(:,jel)) = sp_col.shape_function_gradients(2,1,:,jel);
                end
        end
    end
  end
end
