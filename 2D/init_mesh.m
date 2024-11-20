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
MESH.verts_tri = accumarray(MESH.tri_verts(:), tri_ind(:), [size(MESH.verts, 2), 1], @(A) {A}).';

% edges = [MESH.tri_verts(:,1), MESH.tri_verts(:,2);
%          MESH.tri_verts(:,1), MESH.tri_verts(:,3);
%          MESH.tri_verts(:,2), MESH.tri_verts(:,3)];
% edges = sort(edges, 2);
% edges_adjtri_list = repmat((1:size(MESH.tri_verts, 1)).', [3, 1]);
% [edges, ~, ic] = unique(edges, 'rows');

% MESH.edges = edges;
% MESH.tri_edges = sort(reshape(ic, size(MESH.tri_verts)), 2);

% MESH.verts_valency = accumarray(MESH.edges(:), ones(size(MESH.edges(:)))).';

% edges_adjtri = cell(size(edges, 1), 1);
% for ii = 1:size(edges_adjtri_list, 1)
%     edges_adjtri{ic(ii)}(end+1) = edges_adjtri_list(ii);
% end
% edges_adjtri = cellfun(@(A) resize(A, 2, 'FillValue', nan), edges_adjtri, 'UniformOutput', false);
% edges_adjtri = cell2mat(edges_adjtri);

% MESH.edges_adjtri = edges_adjtri;

end

