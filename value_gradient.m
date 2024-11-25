function [obj,grad] = value_gradient(MESH, opt)

if nargout == 1
    tri_area = tri_area_vec(MESH.tri_verts, MESH.verts);

    if opt.objective_area_squared
        obj = sum(tri_area.^2);
    else
        obj = sum(tri_area);
    end
else
    [tri_area, grad_mat] = tri_area_gradient_vec(MESH.tri_verts, MESH.verts);
    
    if opt.objective_area_squared
        tri_area_temp = shiftdim(tri_area, -1);
        grad_mat = 2 * tri_area_temp .* grad_mat;
    end
    
    dim_ind = repmat((1:3).', [1, 3, size(MESH.tri_verts, 1)]);
    tri_ind = repmat(shiftdim(MESH.tri_verts.', -1), [3, 1, 1]);
    
    grad = accumarray([dim_ind(:), tri_ind(:)], grad_mat(:), [3, size(MESH.verts, 2)]);
    
    if opt.objective_area_squared
        obj = sum(tri_area.^2);
    else
        obj = sum(tri_area);
    end
end


end

