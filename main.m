% clc;
clear;

boundary = [0, 0, 1, 1, 1, 1, 0, 0;
            0, 0, 0, 0, 1, 1, 1, 1;
            0, 1, 1, 0, 0, 1, 1, 0];

% t = linspace(0, 2*pi, 257);
% t(end) = [];
% boundary = [4 * cos(t);
%             4 * sin(t);
%             cos(4 * t)];

% rng(2);
rng('shuffle');
opt.num_refine = 4;
opt.gradient_threshold = 1e-6;
opt.max_iterations = 1e4;
% opt.quality_threshold = @(x) max(0.9 - 1e-4 * x, 0.5);
opt.quality_threshold = 0.75;
opt.subdivide_boundary = true;
opt.verbose = 2;
opt.use_ADAM = true;
opt.gamma = 1e-3;
opt.objective_area_squared = false;
opt.max_swaps = 10;
opt.max_valency = 9;
opt.warmup_iterations = 0;
opt.use_linesearch = false;
[verts, tri] = minsurf3D(boundary, opt);