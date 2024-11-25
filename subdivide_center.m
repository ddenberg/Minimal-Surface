function MESH_new = subdivide_center(MESH)
% The centroid of each triangle is computed and then added to 'verts'.
% Triangles are created by connecting old vertices to the new centers and
% new centers to each other akin to a voronoi tesselation.

MESH_new.verts = [MESH.verts, nan(3, size(MESH.tri_verts,1))];
tri_center_ind = size(MESH.verts, 2)+1:(size(MESH.verts,2)+size(MESH.tri_verts,1));

MESH_new.boundary_verts = false(1, size(MESH_new.verts, 2));
MESH_new.boundary_verts(1:size(MESH.verts,2)) = MESH.boundary_verts;

for ii = 1:size(MESH.tri_verts, 1)
    verts_ii = MESH.verts(:,MESH.tri_verts(ii,:));
    verts_c = mean(verts_ii, 2);
    
    MESH_new.verts(:,size(MESH.verts,2)+ii) = verts_c;
end

old_edges = [MESH.tri_verts(:,1), MESH.tri_verts(:,2);
             MESH.tri_verts(:,1), MESH.tri_verts(:,3);
             MESH.tri_verts(:,2), MESH.tri_verts(:,3)];
old_edges = sort(old_edges, 2);
old_edges = unique(old_edges, 'rows');
old_boundary_edges = all([MESH.boundary_verts(old_edges(:,1)); MESH.boundary_verts(old_edges(:,2))], 1).';

old_edge_adjtri = arrayfun(@(v1, v2) intersect(MESH.verts_tri{v1}, MESH.verts_tri{v2}), ...
    old_edges(:,1), old_edges(:,2), 'UniformOutput', false);

% old_tri_adjtri = cell(size(MESH.tri_verts, 1), 1);
old_tri_edge_ind = cell(size(MESH.tri_verts, 1), 1);
for ii = 1:length(old_edge_adjtri)
    % if length(old_edge_adjtri{ii}) > 1
    %     old_tri_adjtri{old_edge_adjtri{ii}(1)} = union(old_tri_adjtri{old_edge_adjtri{ii}(1)}, old_edge_adjtri{ii}(2));
    %     old_tri_adjtri{old_edge_adjtri{ii}(2)} = union(old_tri_adjtri{old_edge_adjtri{ii}(2)}, old_edge_adjtri{ii}(1));
    % end
    for jj = 1:length(old_edge_adjtri{ii})
        old_tri_edge_ind{old_edge_adjtri{ii}(jj)} = union(old_tri_edge_ind{old_edge_adjtri{ii}(jj)}, ii);
    end
end
old_tri_edge_ind = cell2mat(old_tri_edge_ind);

new_tri_verts = [];
for ii = 1:size(MESH.tri_verts, 1)
    for jj = 1:3 % edges
        if old_boundary_edges(old_tri_edge_ind(ii,jj))
            edge_verts = old_edges(old_tri_edge_ind(ii,jj),:);
            new_tri_verts = [new_tri_verts; sort([edge_verts(1), edge_verts(2), tri_center_ind(ii)])];
        else
            edge_verts = old_edges(old_tri_edge_ind(ii,jj),:);
            neighbor_tri = setdiff(old_edge_adjtri{old_tri_edge_ind(ii,jj)}, ii);
            new_tri_verts = [new_tri_verts; sort([edge_verts(1), tri_center_ind(neighbor_tri), tri_center_ind(ii)])];
            new_tri_verts = [new_tri_verts; sort([edge_verts(2), tri_center_ind(neighbor_tri), tri_center_ind(ii)])];
        end
    end
end
new_tri_verts = unique(new_tri_verts, 'rows');
MESH_new.tri_verts = new_tri_verts;

tri_ind = repmat((1:size(MESH_new.tri_verts, 1)).', [1, 3]);
verts_tri = accumarray(MESH_new.tri_verts(:), tri_ind(:), [size(MESH_new.verts, 2), 1], @(A) {A}).';
MESH_new.verts_tri = cellfun(@(A) sort(A), verts_tri, 'UniformOutput', false);

% MESH_new.tri_timer = randi(250, size(MESH_new.tri_verts, 1), 1);
MESH_new.tri_timer = inf(size(MESH_new.tri_verts, 1), 1);

end

