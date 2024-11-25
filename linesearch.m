function [MESH_new, alpha] = linesearch(p, grad, MESH, opt, obj_old)

MESH_new = MESH;

tau = 0.5;
c = 0.1;

m = grad(:).' * p(:);
t = -c * m;

alpha = 1;

MESH_new.verts = MESH.verts + alpha * p;
obj_new = value_gradient(MESH_new, opt);

while obj_old - obj_new < alpha * t
    alpha = alpha * tau;

    MESH_new.verts = MESH.verts + alpha * p;
    obj_new = value_gradient(MESH_new, opt);
end

end

