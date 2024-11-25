function [A, P, tri_edge_length] = tri_area_perimeter_vec(tri, verts)
% computes the area of a 2D triangle and the gradient of the area with respect
% to each vertex
% INPUT:
%   tri: [m x 3] list of triangles
%   verts: [3 x n] list of coordinates

edge1_vec = verts(:,tri(:,1)) - verts(:,tri(:,2));
edge2_vec = verts(:,tri(:,1)) - verts(:,tri(:,3));
edge3_vec = verts(:,tri(:,2)) - verts(:,tri(:,3));

edge1_length = sqrt(sum(edge1_vec.^2, 1));
edge2_length = sqrt(sum(edge2_vec.^2, 1));
edge3_length = sqrt(sum(edge3_vec.^2, 1));

tri_edge_length = [edge1_length; edge2_length; edge3_length];

P = edge1_length + edge2_length + edge3_length;

detTxy = verts(1,tri(:,1)).*verts(2,tri(:,2)) - verts(1,tri(:,2)).*verts(2,tri(:,1)) - ...
         verts(1,tri(:,1)).*verts(2,tri(:,3)) + verts(1,tri(:,3)).*verts(2,tri(:,1)) + ...
         verts(1,tri(:,2)).*verts(2,tri(:,3)) - verts(1,tri(:,3)).*verts(2,tri(:,2));

detTyz = verts(2,tri(:,1)).*verts(3,tri(:,2)) - verts(2,tri(:,2)).*verts(3,tri(:,1)) - ...
         verts(2,tri(:,1)).*verts(3,tri(:,3)) + verts(2,tri(:,3)).*verts(3,tri(:,1)) + ...
         verts(2,tri(:,2)).*verts(3,tri(:,3)) - verts(2,tri(:,3)).*verts(3,tri(:,2));

detTzx = verts(1,tri(:,2)).*verts(3,tri(:,1)) - verts(1,tri(:,1)).*verts(3,tri(:,2)) + ...
         verts(1,tri(:,1)).*verts(3,tri(:,3)) - verts(1,tri(:,3)).*verts(3,tri(:,1)) - ...
         verts(1,tri(:,2)).*verts(3,tri(:,3)) + verts(1,tri(:,3)).*verts(3,tri(:,2));

A = 1/2 * sqrt(detTxy.^2 + detTyz.^2 + detTzx.^2);

end

