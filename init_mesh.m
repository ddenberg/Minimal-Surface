function MESH = init_mesh(boundary)

MESH.verts = [boundary, mean(boundary, 2)];
MESH.boundary_verts = true(1, size(MESH.verts, 2));
MESH.boundary_verts(size(boundary, 2)+1:end) = false;

tri_verts = nan(size(boundary, 2), 3);
for ii = 1:size(tri_verts, 1)
    tri_verts(ii,:) = [ii, mod(ii,size(tri_verts,1))+1, size(boundary,2)+1];
end
MESH.tri_verts = sort(tri_verts, 2);

tri_ind = repmat((1:size(MESH.tri_verts, 1)).', [1, 3]);
verts_tri = accumarray(MESH.tri_verts(:), tri_ind(:), [size(MESH.verts, 2), 1], @(A) {A}).';
MESH.verts_tri = cellfun(@(A) sort(A), verts_tri, 'UniformOutput', false);

% MESH.tri_timer = randi(250, size(MESH.tri_verts, 1), 1);
MESH.tri_timer = inf(size(MESH.tri_verts, 1), 1);

end

