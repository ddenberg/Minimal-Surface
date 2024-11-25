function MESH_new = subdivide_midpoint(MESH)
% This function subdivides each triangle into 4 new triangles.
% The midpoints of each edge are computed and then new triangles are
% created by connecting the midpoints

midpoint_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

verts_new = MESH.verts;
function midpoint = get_midpoint(v1, v2)

    sorted_edge = sort([v1, v2]);
    key = sprintf('%d-%d', sorted_edge(1), sorted_edge(2));
    
    % If the midpoint is not already computed, calculate it
    if ~isKey(midpoint_map, key)
        midpoint = (verts_new(:,v1) + verts_new(:,v2)) / 2;
        midpoint_map(key) = size(verts_new, 2) + 1;
        verts_new(:,end+1) = midpoint;
    end
    
    midpoint = midpoint_map(key);
end

tri_verts_new = nan(4 * size(MESH.tri_verts, 1), 3);
for ii = 1:size(MESH.tri_verts, 1)
    v1 = MESH.tri_verts(ii, 1);
    v2 = MESH.tri_verts(ii, 2);
    v3 = MESH.tri_verts(ii, 3);
    
    m12 = get_midpoint(v1, v2);
    m23 = get_midpoint(v2, v3);
    m31 = get_midpoint(v3, v1);
    
    % Create the 4 new smaller triangles
    tri_verts_new((ii-1)*4 + 1,:) = [v1, m12, m31];
    tri_verts_new((ii-1)*4 + 2,:) = [v2, m23, m12];
    tri_verts_new((ii-1)*4 + 3,:) = [v3, m31, m23];
    tri_verts_new((ii-1)*4 + 4,:) = [m12, m23, m31];
end

MESH_new.verts = verts_new;
MESH_new.tri_verts = sort(tri_verts_new, 2);

tri_ind = repmat((1:size(MESH_new.tri_verts, 1)).', [1, 3]);
verts_tri = accumarray(MESH_new.tri_verts(:), tri_ind(:), [size(MESH_new.verts, 2), 1], @(A) {A}).';
MESH_new.verts_tri = cellfun(@(A) sort(A), verts_tri, 'UniformOutput', false);

% To find boundary verts look for edge pairs which contain only one
% adjacent triangle

edges = [MESH_new.tri_verts(:,1), MESH_new.tri_verts(:,2);
         MESH_new.tri_verts(:,1), MESH_new.tri_verts(:,3);
         MESH_new.tri_verts(:,2), MESH_new.tri_verts(:,3)];
edges = sort(edges, 2);
edges = unique(edges, 'rows');

edge_adjtri = arrayfun(@(v1, v2) intersect(verts_tri{v1}, verts_tri{v2}), edges(:,1), edges(:,2), 'UniformOutput', false);

% Check that there are no more than 2 and no less than 1 adjacent triangles
% per edge
edge_adjtri_num = cellfun('length', edge_adjtri);
if any(edge_adjtri_num > 2) || any(edge_adjtri_num < 1)
    error('Invalid connectivity found.');
end

boundary_edges = edge_adjtri_num == 1;
boundary_verts = accumarray(edges(:), repmat(boundary_edges(:), [2, 1]), [size(MESH_new.verts, 2), 1]);
MESH_new.boundary_verts = boundary_verts.' > 0;

% MESH_new.tri_timer = randi(250, size(MESH_new.tri_verts, 1), 1);
MESH_new.tri_timer = inf(size(MESH_new.tri_verts, 1), 1);

end

