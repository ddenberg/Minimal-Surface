function MESH_new = subdivide_center_3D(MESH)
% This function subdivides each triangle into 3 new triangles.
% The centroid of each triangle is computed and then added to 'verts'. The
% existing vertices of each triangle are connected to the centroid to form
% 3 new triangles. The old triangle is removed.

MESH_new.verts = [MESH.verts, nan(3, size(MESH.tri_verts,1))];
index_new = size(MESH.verts,2) + 1;

MESH_new.boundary_verts = false(1, size(MESH_new.verts, 2));
MESH_new.boundary_verts(1:size(MESH.verts,2)) = MESH.boundary_verts;

MESH_new.tri_verts = nan(3 * size(MESH.tri_verts,1), 3);
for ii = 1:size(MESH.tri_verts, 1)
    verts_ii = MESH.verts(:,MESH.tri_verts(ii,:));
    verts_c = mean(verts_ii, 2);
    
    MESH_new.verts(:,index_new) = verts_c;
    MESH_new.tri_verts((ii-1)*3 + 1,:) = [MESH.tri_verts(ii,1), MESH.tri_verts(ii,2), index_new];
    MESH_new.tri_verts((ii-1)*3 + 2,:) = [MESH.tri_verts(ii,1), MESH.tri_verts(ii,3), index_new];
    MESH_new.tri_verts((ii-1)*3 + 3,:) = [MESH.tri_verts(ii,2), MESH.tri_verts(ii,3), index_new];
    index_new = index_new + 1;
end
MESH_new.tri_verts = sort(MESH_new.tri_verts, 2);

tri_ind = repmat((1:size(MESH_new.tri_verts, 1)).', [1, 3]);
MESH_new.verts_tri = accumarray(MESH_new.tri_verts(:), tri_ind(:), [size(MESH_new.verts, 2), 1], @(A) {A}).';

end

