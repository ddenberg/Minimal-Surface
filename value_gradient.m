function [obj, grad] = value_gradient(MESH, opt)

% tri_area = zeros(size(MESH.tri_verts, 1), 1);
% % tri_perim = zeros(size(MESH.tri_verts, 1), 1);
% grad = zeros(size(MESH.verts));
% for ii = 1:size(MESH.tri_verts, 1)
%     verts_ii = MESH.verts(:,MESH.tri_verts(ii,:));
%     [A, gradA] = tri_area_gradient_3D(verts_ii);
% 
%     tri_area(ii) = A;
%     % tri_perim(ii) = P;
% 
%     % NOTE: This can potentially be done faster by writing gradA to a [size(triangles,1), 3] array
%     % and then using 'accumarray'
%     if opt.objective_area_squared
%         grad(:,MESH.tri_verts(ii,:)) = grad(:,MESH.tri_verts(ii,:)) + 2 * A * gradA;
%     else
%         grad(:,MESH.tri_verts(ii,:)) = grad(:,MESH.tri_verts(ii,:)) + gradA;
%     end
% 
% end

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

