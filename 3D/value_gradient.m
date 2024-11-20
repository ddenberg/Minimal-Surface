function [obj, grad] = value_gradient(MESH, opt)

tri_areas = zeros(size(MESH.tri_verts, 1), 1);
grad = zeros(size(MESH.verts));
for ii = 1:size(MESH.tri_verts, 1)
    verts_ii = MESH.verts(:,MESH.tri_verts(ii,:));
    [A, gradA] = tri_area_gradient_3D(verts_ii);

    tri_areas(ii) = A;

    % NOTE: This can potentially be done faster by writing gradA to a [size(triangles,1), 3] array
    % and then using 'accumarray'
    if opt.objective_area_squared
        grad(:,MESH.tri_verts(ii,:)) = grad(:,MESH.tri_verts(ii,:)) + 2 * A * gradA;
    else
        grad(:,MESH.tri_verts(ii,:)) = grad(:,MESH.tri_verts(ii,:)) + gradA;
    end

end

if opt.objective_area_squared
    obj = sum(tri_areas.^2);
else
    obj = sum(tri_areas);
end


end

