clc;
clear;

% tri = [1, 2, 3];
x = sym('x', [3, 3]);
assume(x, 'real');

Txy = [x([1,2],:); ones(1,3)];
Tyz = [x([2,3],:); ones(1,3)];
Tzx = [x([3,1],:); ones(1,3)];

detTxy = det(Txy);
detTyz = det(Tyz);
detTzx = det(Tzx);

A = 1/2 * sqrt(detTxy^2 + detTyz^2 + detTzx^2);

dA_ddetTxy = detTxy / (4 * A);
dA_ddetTyz = detTyz / (4 * A);
dA_ddetTzx = detTzx / (4 * A);

ddetTxy_dx = reshape(jacobian(detTxy, x(:)), size(x));
ddetTyz_dx = reshape(jacobian(detTyz, x(:)), size(x));
ddetTzx_dx = reshape(jacobian(detTzx, x(:)), size(x));

dA_dx = dA_ddetTxy * ddetTxy_dx + dA_ddetTyz * ddetTyz_dx + dA_ddetTzx * ddetTzx_dx;

% ddetTxy_dx = [x2_2 - x2_3, x2_3 - x2_1, x2_1 - x2_2]
%              [x1_3 - x1_2, x1_1 - x1_3, x1_2 - x1_1]
%              [          0,           0,           0]

% ddetTyz_dx = [          0,           0,           0]
%              [x3_2 - x3_3, x3_3 - x3_1, x3_1 - x3_2]
%              [x2_3 - x2_2, x2_1 - x2_3, x2_2 - x2_1]

% ddetTzx_dx = [x3_3 - x3_2, x3_1 - x3_3, x3_2 - x3_1]
%              [          0,           0,           0]
%              [x1_2 - x1_3, x1_3 - x1_1, x1_1 - x1_2]