clc;
clear;

% boundary = [0, 0, 1, 1, 1, 1, 0, 0;
%             0, 0, 0, 0, 1, 1, 1, 1;
%             0, 1, 1, 0, 0, 1, 1, 0];

t = linspace(0, 2*pi, 101);
t(end) = [];
boundary = [4 * cos(t);
            4 * sin(t);
            cos(4 * t)];

rng(2);
opt.num_refine = 3;
opt.gradient_threshold = 1e-7;
opt.max_iterations = 1e4;
% opt.quality_threshold = @(x) max(0.9 - 1e-4 * x, 0.5);
opt.quality_threshold = 0.8;
opt.subdivide_boundary = false;
opt.verbose = 2;
opt.use_ADAM = false;
opt.gamma = 1e-2;
opt.objective_area_squared = true;
opt.max_swaps = 5;
opt.max_valency = 8;
[verts, tri] = minsurf3D(boundary, opt);

