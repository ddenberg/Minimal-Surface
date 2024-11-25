function MESH_new = remesh(MESH)

MESH_new = MESH;

tri_verts = delaunay(MESH_new.verts(1,:), MESH_new.verts(2,:));
MESH_new.tri_verts = tri_verts;

tri_ind = repmat((1:size(MESH_new.tri_verts, 1)).', [1, 3]);
verts_tri = accumarray(MESH_new.tri_verts(:), tri_ind(:), [size(MESH_new.verts, 2), 1], @(A) {A}).';
MESH_new.verts_tri = cellfun(@(A) sort(A), verts_tri, 'UniformOutput', false);

end

