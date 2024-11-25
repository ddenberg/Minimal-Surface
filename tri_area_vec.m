function A = tri_area_vec(tri, verts)
% computes the area of a 3D triangle
% INPUT:
%   tri: [m x 3] list of triangles
%   verts: [3 x n] list of coordinates

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

