function MESH_new = subdivide_center_2D(MESH)
% This function subdivides each triangle into 3 new triangles.
% The centroid of each triangle is computed and then added to 'verts'. The
% existing vertices of each triangle are connected to the centroid to form
% 3 new triangles. The old triangle is removed.

MESH_new.verts = [MESH.verts, nan(2, size(MESH.tri_verts,1))];
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

% edges = [MESH_new.tri_verts(:,1), MESH_new.tri_verts(:,2);
%          MESH_new.tri_verts(:,1), MESH_new.tri_verts(:,3);
%          MESH_new.tri_verts(:,2), MESH_new.tri_verts(:,3)];
% edges = sort(edges, 2);
% edges_adjtri_list_new = repmat((1:size(MESH_new.tri_verts, 1)).', [3, 1]);
% [edges, ~, ic] = unique(edges, 'rows');
% 
% MESH_new.edges = edges;
% MESH_new.tri_edges = sort(reshape(ic, size(MESH_new.tri_verts)), 2);
% 
% MESH_new.verts_valency = accumarray(MESH_new.edges(:), ones(size(MESH_new.edges(:)))).';
% 
% edges_adjtri = cell(size(edges, 1), 1);
% for ii = 1:size(edges_adjtri_list_new, 1)
%     edges_adjtri{ic(ii)}(end+1) = edges_adjtri_list_new(ii);
% end
% edges_adjtri = cellfun(@(A) resize(A, 2, 'FillValue', nan), edges_adjtri, 'UniformOutput', false);
% edges_adjtri = cell2mat(edges_adjtri);
% 
% MESH_new.edges_adjtri = edges_adjtri;

end

