clc;
clear;

tri = [1, 2, 3];
x = sym('x', [2, 3]);
assume(x, 'real');

T = [x(1,tri(1,:));
     x(2,tri(1,:));
     ones(1,3)];

detT = det(T);
A = 1/2 * abs(detT);
Asq = (1/2 * detT).^2;

% A = abs(x1_1*x2_2 - x1_2*x2_1 - x1_1*x2_3 + x1_3*x2_1 + x1_2*x2_3 - x1_3*x2_2)/2

dA_dx = reshape(jacobian(A, x(:)), size(x));
dAsq_dx = reshape(jacobian(Asq, x(:)), size(x));

% ddetT_dx = [x2_2 - x2_3, x2_3 - x2_1, x2_1 - x2_2;
%             x1_3 - x1_2, x1_1 - x1_3, x1_2 - x1_1];

% dA_dx = 1/2*sign(detT) * ddetT_dx