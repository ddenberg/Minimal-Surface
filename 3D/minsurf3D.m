function [verts, tri] = minsurf3D(boundary, opt)
% minsurf computes the minimal surface defined by a closed boundary
% INPUTS:
%   boundary: [3 x n] ordered list of points on the boundary curve
%   opt.num_refine: number of triangle subdivisions to use
%   opt.gradient_threshold: controls when to stop gradient descent
%   opt.max_iterations: maximum number of iterations per refinement
%   opt.quality_threshold: triangle quality threshold. Closer to 1 is more
%       equilateral.
%   opt.subdivide_boundary: (true/false) Set to true if you are providing a
%       coarse boundary and false if your boundary is already highly
%       subdivided.
%   opt.use_ADAM: (true/false) Set to true to use ADAM optimizer and false
%       for regular gradient descent. (NOTE: YOU MAY ENCOUNTER JITTERING
%       WITH ADAM)
%   opt.gamma: the step size for ADAM or gradient descent. If things blow
%       up, consider first decreasing gamma.
%   opt.objective_area_squared: Use sum(tri_area^2) rather than
%       sum(tri_area) for the objective.
%   opt.verbose: true/false flag for displaying progress

MESH = init_mesh(boundary);

% ADAM parameters
beta1 = 0.9;
beta2 = 0.999;
gamma = opt.gamma;
epsilon = 1e-8;

for refine = 0:opt.num_refine
    if refine > 0
        if opt.subdivide_boundary
            MESH = subdivide_midpoint_3D(MESH);
            allow_boundary_swaps = false;
        else
            MESH = subdivide_center_3D(MESH);
            allow_boundary_swaps = true;
        end
        
    else
        allow_boundary_swaps = false;
    end

    MESH = refine_triangles(MESH, opt.quality_threshold, allow_boundary_swaps);

    [obj, grad] = value_gradient(MESH, opt);
    m = zeros(size(MESH.verts));
    v = zeros(size(MESH.verts));
    
    
    grad_norm = rms(grad(:,~MESH.boundary_verts), 'all');
    grad_norm_tracker = grad_norm;
    obj_tracker = obj;
    iter = 1;
    while grad_norm > opt.gradient_threshold && iter < opt.max_iterations
        if opt.use_ADAM
            m = beta1 * m + (1 - beta1) * grad;
            v = beta2 * v + (1 - beta2) * grad.^2;
            m_hat = m / (1 - beta1^iter);
            v_hat = v / (1 - beta2^iter);
    
            delta = -gamma * m_hat ./ (sqrt(v_hat) + epsilon);
            MESH.verts(:,~MESH.boundary_verts) = MESH.verts(:,~MESH.boundary_verts) + ...
                delta(:,~MESH.boundary_verts);
            m(:,MESH.boundary_verts) = 0;
            v(:,MESH.boundary_verts) = 0;
        else
            MESH.verts(:,~MESH.boundary_verts) = MESH.verts(:,~MESH.boundary_verts) - ...
                gamma * grad(:,~MESH.boundary_verts);
        end

        if mod(iter, 1) == 0
            MESH = refine_triangles(MESH, opt.quality_threshold, allow_boundary_swaps);
        end
        
        [obj, grad] = value_gradient(MESH, opt);

        grad_norm = rms(grad(:,~MESH.boundary_verts), 'all');
        grad_norm_tracker(end+1) = grad_norm;
        obj_tracker(end+1) = obj;

        if mod(iter, 100) == 0 && opt.verbose > 1
            figure(1);
            cla;
            hold on;
            trimesh(MESH.tri_verts, MESH.verts(1,:), MESH.verts(2,:), MESH.verts(3,:), 'EdgeColor', 'k');
            scatter3(MESH.verts(1,:), MESH.verts(2,:), MESH.verts(3,:), 'filled');
            axis equal vis3d;
            xlim([min(boundary(1,:)) - 0.25, max(boundary(1,:)) + 0.25]);
            ylim([min(boundary(2,:)) - 0.25, max(boundary(2,:)) + 0.25]); 
            zlim([min(boundary(3,:)) - 0.25, max(boundary(3,:)) + 0.25]);
            title(sprintf('iter = %d, grad norm = %e', iter, grad_norm));
            drawnow;
        end

        iter = iter + 1;
    end

    if opt.verbose > 0
        fprintf('Done refine level: %d, Num triangles: %d, gradient norm: %e, iterations: %d\n', ...
            refine, size(MESH.tri_verts, 1), grad_norm, iter);
    end

    if opt.verbose > 1
        figure(1);
        cla;
        hold on;
        trimesh(MESH.tri_verts, MESH.verts(1,:), MESH.verts(2,:), MESH.verts(3,:), 'EdgeColor', 'k');
        scatter3(MESH.verts(1,:), MESH.verts(2,:), MESH.verts(3,:), 'filled');
        axis equal vis3d;
        xlim([min(boundary(1,:)) - 0.25, max(boundary(1,:)) + 0.25]);
        ylim([min(boundary(2,:)) - 0.25, max(boundary(2,:)) + 0.25]); 
        zlim([min(boundary(3,:)) - 0.25, max(boundary(3,:)) + 0.25]);  
        title(sprintf('iter = %d, grad norm = %e', iter, grad_norm));
        drawnow;
    end

    if opt.verbose > 1
        figure(2);
        subplot(opt.num_refine+1, 2, refine*2 + 1);
        hold on;
        plot(obj_tracker, 'LineWidth', 2);
        set(gca, 'xscale', 'log');
        ylabel('Objective');

        subplot(opt.num_refine+1, 2, refine*2 + 2);
        hold on;
        plot(grad_norm_tracker, 'LineWidth', 2);
        set(gca, 'xscale', 'log', 'yscale', 'log');
        ylabel('Gradient Norm');
    end
end

verts = MESH.verts;
tri = MESH.tri_verts;


end
