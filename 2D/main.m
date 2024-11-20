clc;
clear;

boundary = [0, 1, 1.5, 1, 0;
            0, 0, 0.5, 1, 1];
% boundary = [cos(linspace(0, pi, 20)), linspace(-1, 1, 10)
%             sin(linspace(0, pi, 20)), zeros(1, 10)];

opt.num_refine = 2;
opt.gradient_threshold = 1e-8;
opt.max_iterations = 1e4;
opt.quality_threshold = 0.6;
opt.subdivide_boundary = true;
opt.verbose = 2;
[verts, tri] = minsurf2D(boundary, opt);

