function [A, gradA] = tri_area_gradient_2D(verts)
% computes the area of a 2D triangle and the gradient of the area with respect
% to each vertex
% INPUT:
%   x: [1 x 3] list of x coordinates
%   y: [1 x 3] list of y coordinates

T = [verts; ones(1,3)];
detT = det(T);
A = 1/2 * abs(detT);

gradT = [verts(2,2) - verts(2,3), verts(2,3) - verts(2,1), verts(2,1) - verts(2,2);
         verts(1,3) - verts(1,2), verts(1,1) - verts(1,3), verts(1,2) - verts(1,1)];

gradA = 1/2 * sign(detT) * gradT;

end

