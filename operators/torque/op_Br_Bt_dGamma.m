% op_Br_Bt_dGamma: calculates the torque resulting from a b_field
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points


function T = op_Br_Bt_dGamma (space, msh, Az)

    %shpu_derx = reshape (space.shape_function_gradients(1,:,:,:), space.ncomp, msh.nqn, space.nsh_max, msh.nel);
    %shpu_dery = reshape (space.shape_function_gradients(2,:,:,:), space.ncomp, msh.nqn, space.nsh_max, msh.nel);

    gradu = reshape (space.shape_function_gradients, space.ncomp, [], msh.nqn, space.nsh_max, msh.nel);
    
    jacdet_weights = reshape (msh.jacdet .* msh.quad_weights, [1, msh.nqn, msh.nel]);

    ndim = size (gradu, 2);
    
    T = 0;
    
    for iel = 1:msh.nel
        if (all (msh.jacdet(:,iel)))
            %shpu_derx_iel = reshape (shpu_derx(:, :, :, iel), space.ncomp, msh.nqn, space.nsh_max, 1);
            %shpu_dery_iel = reshape (shpu_dery(:, :, :, iel), space.ncomp, msh.nqn, space.nsh_max, 1);

            gradu_iel = gradu(:,:,:,:,iel);
            gradu_iel_x = reshape(gradu_iel(:,1,:,:), space.ncomp, msh.nqn, space.nsh_max);
            gradu_iel_y = reshape(gradu_iel(:,2,:,:), space.ncomp, msh.nqn, space.nsh_max);

            Az_iel = reshape(Az(space.connectivity(:,iel)), [1, 1, space.nsh_max]);
            normal_iel = reshape (msh.normal(:,:,iel), [ndim, msh.nqn]);

            Ni_x_uiel = sum(gradu_iel_x.*Az_iel, 3);
            Ni_y_uiel = sum(gradu_iel_y.*Az_iel, 3);
            
            bx = Ni_y_uiel;
            by = -Ni_x_uiel;
            Br = bx.*normal_iel(1,:) + by.*normal_iel(2,:);
            Bt = -bx.*normal_iel(2,:) + by.*normal_iel(1,:);
            r = reshape((msh.geo_map(1,:,iel).^2 + msh.geo_map(2,:,iel).^2).^0.5, [1, msh.nqn]);
%             bmag = (bx.^2 + by.^2).^0.5;
%             bmag1 = (Br.^2 + Bt.^2).^0.5;
            %quiver(msh.geo_map(1,:,iel), msh.geo_map(2,:,iel), bx, by);

            jacdet_iel = reshape (jacdet_weights(:,:,iel), [1,msh.nqn]);

            T = T + sum(r.*Br.*Bt.*jacdet_iel, 2);
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_Br_Bt_dGamma: singular map in element number %d', iel)
        end
    end

end


