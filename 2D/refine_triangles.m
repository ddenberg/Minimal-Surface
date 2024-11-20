function MESH_new = refine_triangles(MESH, quality_threshold, allow_boundary_swaps)

MESH_new = MESH;

for ii = 1:size(MESH_new.tri_verts, 1)
    verts_ind = MESH_new.tri_verts(ii,:);

    edges_ii = [verts_ind(1), verts_ind(2);
                verts_ind(1), verts_ind(3);
                verts_ind(2), verts_ind(3)];

    edges_ii_length = vecnorm(MESH_new.verts(:,edges_ii(:,1)) - ...
                              MESH_new.verts(:,edges_ii(:,2)), 2, 1);
    % edges_ii_length = vecnorm(MESH_new.verts(:,MESH_new.edges(edges_ind,1)) - ...
    %                           MESH_new.verts(:,MESH_new.edges(edges_ind,2)), 2, 1);

    T = [MESH_new.verts(:,verts_ind); ones(1, 3)];
    area = 0.5 * abs(det(T));
    perimeter = sum(edges_ii_length);

    quality = 36 * area / (sqrt(3) * perimeter^2);

    if quality < quality_threshold
        % fprintf('Found bad triangle\n');

        [~, longest_edge_ind] = max(edges_ii_length);
        longest_edge = edges_ii(longest_edge_ind,:);

        longest_edge_verts_adjtri = MESH_new.verts_tri(longest_edge);
        longest_edge_adjtri = intersect(longest_edge_verts_adjtri{:});

        if length(longest_edge_adjtri) < 2
            % Can't swap edge in this case because the edge is not adjacent
            % to two triangles.
            continue;
        end

        new_edge = setdiff(MESH_new.tri_verts(longest_edge_adjtri,:), longest_edge);
        new_edge = new_edge(:).';

        longest_edge_verts_valency = cellfun('length', longest_edge_verts_adjtri);
        filter_valency = longest_edge_verts_valency > 4;
        filter_boundary = MESH_new.boundary_verts(longest_edge);
        filter_proceed = filter_valency | filter_boundary;
        if ~all(filter_proceed)
            % If the valency of either of the vertices in the longest edge
            % is less than 4 then don't remove an edge. Unless the vertex
            % with valency less than 4 is a boundary vertex.
            continue;
        end

        if all(MESH_new.boundary_verts(new_edge)) && ~allow_boundary_swaps
            % If both opposite vertices are boundary vertices, then don't
            % do the swap. In the future maybe remove this restriction. The
            % reason would be to avoid very small triangles (poor quality)
            % at the boundary if the boundary is very subdivided.
            continue;
        end


        MESH_old = MESH_new;
        MESH_new.tri_verts(longest_edge_adjtri(1),:) = sort([new_edge, longest_edge(1)]);
        MESH_new.tri_verts(longest_edge_adjtri(2),:) = sort([new_edge, longest_edge(2)]);

        % MESH_new.verts_tri{longest_edge(1)} = setdiff(MESH_new.verts_tri{longest_edge(1)}, ...
        %     longest_edge_adjtri(2));
        % MESH_new.verts_tri{longest_edge(2)} = setdiff(MESH_new.verts_tri{longest_edge(2)}, ...
        %     longest_edge_adjtri(1));
        % 
        % MESH_new.verts_tri{new_edge(1)} = union(MESH_new.verts_tri{new_edge(1)}, ...
        %     longest_edge_adjtri);
        % MESH_new.verts_tri{new_edge(2)} = union(MESH_new.verts_tri{new_edge(2)}, ...
        %     longest_edge_adjtri);

        tri_ind = repmat((1:size(MESH_new.tri_verts, 1)).', [1, 3]);
        MESH_new.verts_tri = accumarray(MESH_new.tri_verts(:), tri_ind(:), [size(MESH_new.verts, 2), 1], @(A) {A}).';


        % check to make sure that the new edge does not have more than 2
        % and less than 1 adjacent triangles
        new_edge_verts_adjtri = MESH_new.verts_tri(new_edge);
        new_edge_adjtri = intersect(new_edge_verts_adjtri{:});
        if length(new_edge_adjtri) > 2 || length(new_edge_adjtri) < 1
            error('Error in new edge connectivity.');
        end

        % figure(1);
        % cla;
        % hold on;
        % triplot(MESH_new.tri_verts, MESH_new.verts(1,:), MESH_new.verts(2,:), 'Color', 'blue');
        % triplot(MESH_new.tri_verts(adjacent_tri,:), MESH_new.verts(1,:), MESH_new.verts(2,:), 'Color', 'red', 'LineWidth', 4);
        % 
        % % tri 1
        % edges_temp = MESH_new.edges(MESH_new.tri_edges(adjacent_tri(1),:),:);
        % for jj = 1:size(edges_temp, 1)
        %     plot(MESH.verts(1, edges_temp(jj,:)), MESH.verts(2, edges_temp(jj,:)), 'k', 'LineWidth', 1);
        % end
        % 
        % % pause;
        % figure(1);
        % cla;
        % hold on;
        % triplot(MESH_new.tri_verts, MESH_new.verts(1,:), MESH_new.verts(2,:), 'Color', 'blue');
        % triplot(MESH_new.tri_verts(adjacent_tri,:), MESH_new.verts(1,:), MESH_new.verts(2,:), 'Color', 'red', 'LineWidth', 4);

        % % tri 2
        % edges_temp = MESH_new.edges(MESH_new.tri_edges(adjacent_tri(2),:),:);
        % for jj = 1:size(edges_temp, 1)
        %     plot(MESH.verts(1, edges_temp(jj,:)), MESH.verts(2, edges_temp(jj,:)), 'k', 'LineWidth', 1);
        % end
        
        % pause;

        rest_ind = setdiff((1:size(MESH_new.tri_verts, 1)).', longest_edge_adjtri);
        [TF, err_ind] = ismember(MESH_new.tri_verts(longest_edge_adjtri,:), ...
                                 MESH_new.tri_verts(rest_ind,:), 'rows');
        if any(TF)
            error('ruh roh');
        end

        [TF, err_ind] = ismember(MESH_new.tri_verts(longest_edge_adjtri(1),:), ...
                                 MESH_new.tri_verts(longest_edge_adjtri(2),:), 'rows');
        if any(TF)
            error('ruh roh');
        end
        

    end

end

