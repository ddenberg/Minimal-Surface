function [A, gradA] = tri_area_gradient_vec(tri, verts)
% computes the area of a 2D triangle and the gradient of the area with respect
% to each vertex
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

dA_ddetTxy = detTxy ./ (4 * A);
dA_ddetTyz = detTyz ./ (4 * A);
dA_ddetTzx = detTzx ./ (4 * A);

ddetTxy_dx = zeros(3, 3, size(tri, 1));
ddetTxy_dx(1,1,:) = verts(2,tri(:,2)) - verts(2,tri(:,3));
ddetTxy_dx(1,2,:) = verts(2,tri(:,3)) - verts(2,tri(:,1));
ddetTxy_dx(1,3,:) = verts(2,tri(:,1)) - verts(2,tri(:,2));
ddetTxy_dx(2,1,:) = verts(1,tri(:,3)) - verts(1,tri(:,2));
ddetTxy_dx(2,2,:) = verts(1,tri(:,1)) - verts(1,tri(:,3));
ddetTxy_dx(2,3,:) = verts(1,tri(:,2)) - verts(1,tri(:,1));

ddetTyz_dx = zeros(3, 3, size(tri, 1));
ddetTyz_dx(2,1,:) = verts(3,tri(:,2)) - verts(3,tri(:,3));
ddetTyz_dx(2,2,:) = verts(3,tri(:,3)) - verts(3,tri(:,1));
ddetTyz_dx(2,3,:) = verts(3,tri(:,1)) - verts(3,tri(:,2));
ddetTyz_dx(3,1,:) = verts(2,tri(:,3)) - verts(2,tri(:,2));
ddetTyz_dx(3,2,:) = verts(2,tri(:,1)) - verts(2,tri(:,3));
ddetTyz_dx(3,3,:) = verts(2,tri(:,2)) - verts(2,tri(:,1));

ddetTzx_dx = zeros(3, 3, size(tri, 1));
ddetTzx_dx(1,1,:) = verts(3,tri(:,3)) - verts(3,tri(:,2));
ddetTzx_dx(1,2,:) = verts(3,tri(:,1)) - verts(3,tri(:,3));
ddetTzx_dx(1,3,:) = verts(3,tri(:,2)) - verts(3,tri(:,1));
ddetTzx_dx(3,1,:) = verts(1,tri(:,2)) - verts(1,tri(:,3));
ddetTzx_dx(3,2,:) = verts(1,tri(:,3)) - verts(1,tri(:,1));
ddetTzx_dx(3,3,:) = verts(1,tri(:,1)) - verts(1,tri(:,2));

% ddetTxy_dx = [verts(2,2) - verts(2,3), verts(2,3) - verts(2,1), verts(2,1) - verts(2,2);
%               verts(1,3) - verts(1,2), verts(1,1) - verts(1,3), verts(1,2) - verts(1,1);
%                                     0,                       0,                       0];
% 
% ddetTyz_dx = [                      0,                       0,                       0;
%               verts(3,2) - verts(3,3), verts(3,3) - verts(3,1), verts(3,1) - verts(3,2);
%               verts(2,3) - verts(2,2), verts(2,1) - verts(2,3), verts(2,2) - verts(2,1)];
% 
% ddetTzx_dx = [verts(3,3) - verts(3,2), verts(3,1) - verts(3,3), verts(3,2) - verts(3,1);
%                                     0,                       0,                       0;
%               verts(1,2) - verts(1,3), verts(1,3) - verts(1,1), verts(1,1) - verts(1,2)];

dA_ddetTxy = reshape(dA_ddetTxy, [1, 1, size(dA_ddetTxy, 2)]);
dA_ddetTyz = reshape(dA_ddetTyz, [1, 1, size(dA_ddetTyz, 2)]);
dA_ddetTzx = reshape(dA_ddetTzx, [1, 1, size(dA_ddetTzx, 2)]);

gradA = dA_ddetTxy .* ddetTxy_dx + dA_ddetTyz .* ddetTyz_dx + dA_ddetTzx .* ddetTzx_dx;

end

