clc;
clear;

boundary = [0, 0, 1, 1, 1, 1, 0, 0;
            0, 0, 0, 0, 1, 1, 1, 1;
            0, 1, 1, 0, 0, 1, 1, 0];

opt.num_refine = 3;
opt.gradient_threshold = 1e-6;
opt.max_iterations = 1e4;
opt.quality_threshold = 0.75;
opt.subdivide_boundary = true;
opt.verbose = 2;
opt.use_ADAM = false;
opt.gamma = 1e-1;
opt.objective_area_squared = false;
[verts, tri] = minsurf3D(boundary, opt);

