function [A, gradA, P] = tri_area_gradient_3D(verts)
% computes the area of a 2D triangle and the gradient of the area with respect
% to each vertex
% INPUT:
%   verts: [3 x 3] list of coordinates

edge1_vec = verts(:,1) - verts(:,2);
edge2_vec = verts(:,1) - verts(:,3);
edge3_vec = verts(:,2) - verts(:,3);

edge1_length = norm(edge1_vec);
edge2_length = norm(edge2_vec);
edge3_length = norm(edge3_vec);

P = edge1_length + edge2_length + edge3_length;

% Txy = [verts([1,2],:); ones(1,3)];
% Tyz = [verts([2,3],:); ones(1,3)];
% Tzx = [verts([3,1],:); ones(1,3)];

% detTxy = det(Txy);
% detTyz = det(Tyz);
% detTzx = det(Tzx);

detTxy = verts(1,1)*verts(2,2) - verts(1,2)*verts(2,1) - verts(1,1)*verts(2,3) + ...
         verts(1,3)*verts(2,1) + verts(1,2)*verts(2,3) - verts(1,3)*verts(2,2);

detTyz = verts(2,1)*verts(3,2) - verts(2,2)*verts(3,1) - verts(2,1)*verts(3,3) + ...
         verts(2,3)*verts(3,1) + verts(2,2)*verts(3,3) - verts(2,3)*verts(3,2);

detTzx = verts(1,2)*verts(3,1) - verts(1,1)*verts(3,2) + verts(1,1)*verts(3,3) - ...
         verts(1,3)*verts(3,1) - verts(1,2)*verts(3,3) + verts(1,3)*verts(3,2);

A = 1/2 * sqrt(detTxy^2 + detTyz^2 + detTzx^2);

dA_ddetTxy = detTxy / (4 * A);
dA_ddetTyz = detTyz / (4 * A);
dA_ddetTzx = detTzx / (4 * A);

ddetTxy_dx = [verts(2,2) - verts(2,3), verts(2,3) - verts(2,1), verts(2,1) - verts(2,2);
              verts(1,3) - verts(1,2), verts(1,1) - verts(1,3), verts(1,2) - verts(1,1);
                                    0,                       0,                       0];

ddetTyz_dx = [                      0,                       0,                       0;
              verts(3,2) - verts(3,3), verts(3,3) - verts(3,1), verts(3,1) - verts(3,2);
              verts(2,3) - verts(2,2), verts(2,1) - verts(2,3), verts(2,2) - verts(2,1)];

ddetTzx_dx = [verts(3,3) - verts(3,2), verts(3,1) - verts(3,3), verts(3,2) - verts(3,1);
                                    0,                       0,                       0;
              verts(1,2) - verts(1,3), verts(1,3) - verts(1,1), verts(1,1) - verts(1,2)];

gradA = dA_ddetTxy * ddetTxy_dx + dA_ddetTyz * ddetTyz_dx + dA_ddetTzx * ddetTzx_dx;

% ddetTxy_dx = [x2_2 - x2_3, x2_3 - x2_1, x2_1 - x2_2]
%              [x1_3 - x1_2, x1_1 - x1_3, x1_2 - x1_1]
%              [          0,           0,           0]

% ddetTyz_dx = [          0,           0,           0]
%              [x3_2 - x3_3, x3_3 - x3_1, x3_1 - x3_2]
%              [x2_3 - x2_2, x2_1 - x2_3, x2_2 - x2_1]

% ddetTzx_dx = [x3_3 - x3_2, x3_1 - x3_3, x3_2 - x3_1]
%              [          0,           0,           0]
%              [x1_2 - x1_3, x1_3 - x1_1, x1_1 - x1_2]

end

