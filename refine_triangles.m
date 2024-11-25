function MESH_new = refine_triangles(MESH, quality_threshold, max_swaps, max_valency)

MESH_new = MESH;
MESH_new.tri_timer = MESH_new.tri_timer + 1;


% tri_order = randperm(size(MESH_new.tri_verts, 1));
% tri_area = nan(1, size(MESH_new.tri_verts, 1));
% tri_perim = nan(1, size(MESH_new.tri_verts, 1));
% tri_edge_length = nan(3, size(MESH_new.tri_verts, 1));
% tri_quality = 36 * tri_area ./ (sqrt(3) * tri_perim.^2);
% recompute_tri_quality = true(length(tri_quality), 1);

[tri_area, tri_perim, tri_edge_length] = tri_area_perimeter_vec(MESH_new.tri_verts, MESH_new.verts);
tri_quality = 36 * tri_area ./ (sqrt(3) * tri_perim.^2);
recompute_tri_quality = false(length(tri_quality), 1);
[~, tri_order] = sort(tri_quality, 'ascend');

swap_count = 0;
for ii = 1:length(tri_order)
    if swap_count >= max_swaps
        break;
    end

    tri_ind = tri_order(ii);

    % if MESH_new.tri_timer(tri_ind) < 250
    %     continue;
    % end
    verts_ind = MESH_new.tri_verts(tri_ind,:);

    edges_ii = [verts_ind(1), verts_ind(2);
                verts_ind(1), verts_ind(3);
                verts_ind(2), verts_ind(3)];

    verts = MESH_new.verts(:,verts_ind);

    if recompute_tri_quality(tri_ind)
        edges_ii_length = sqrt(sum([(verts(:,1) - verts(:,2)).^2, ...
                                    (verts(:,1) - verts(:,3)).^2, ...
                                    (verts(:,2) - verts(:,3)).^2], 1));

        tri_edge_length(:,tri_ind) = edges_ii_length;
    
        detTxy = verts(1,1)*verts(2,2) - verts(1,2)*verts(2,1) - verts(1,1)*verts(2,3) + ...
                 verts(1,3)*verts(2,1) + verts(1,2)*verts(2,3) - verts(1,3)*verts(2,2);
    
        detTyz = verts(2,1)*verts(3,2) - verts(2,2)*verts(3,1) - verts(2,1)*verts(3,3) + ...
                 verts(2,3)*verts(3,1) + verts(2,2)*verts(3,3) - verts(2,3)*verts(3,2);
        
        detTzx = verts(1,2)*verts(3,1) - verts(1,1)*verts(3,2) + verts(1,1)*verts(3,3) - ...
                 verts(1,3)*verts(3,1) - verts(1,2)*verts(3,3) + verts(1,3)*verts(3,2);
    
        area = 1/2 * sqrt(detTxy^2 + detTyz^2 + detTzx^2);
        perimeter = sum(edges_ii_length);
    
        quality = 36 * area / (sqrt(3) * perimeter^2);

        tri_area(tri_ind) = area;
        tri_perim(tri_ind) = perimeter;
        tri_quality(tri_ind) = quality;

        recompute_tri_quality(tri_ind) = false;
    else
        edges_ii_length = tri_edge_length(:,tri_ind);
        quality = tri_quality(tri_ind);
    end

    if quality < quality_threshold
        % fprintf('Found bad triangle\n');

        [~, longest_edge_ind] = max(edges_ii_length);
        longest_edge = edges_ii(longest_edge_ind,:);

        longest_edge_verts_adjtri = MESH_new.verts_tri(longest_edge);
        % TF = ismember(longest_edge_verts_adjtri{1}, longest_edge_verts_adjtri{2});
        TF = ismembc(longest_edge_verts_adjtri{1}, longest_edge_verts_adjtri{2});
        longest_edge_adjtri = longest_edge_verts_adjtri{1}(TF);
        % longest_edge_adjtri = intersect(longest_edge_verts_adjtri{:});

        if length(longest_edge_adjtri) < 2
            % Can't swap edge in this case because the edge is not adjacent
            % to two triangles.
            continue;
        end

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

        all_verts = MESH_new.tri_verts(longest_edge_adjtri,:);
        all_verts = all_verts(:);

        new_edge = all_verts(~(all_verts == longest_edge(1) | all_verts == longest_edge(2)));
        % new_edge = setdiff(MESH_new.tri_verts(longest_edge_adjtri,:), longest_edge);
        new_edge = new_edge(:).';

        if all(MESH_new.boundary_verts(new_edge))
            % If both opposite vertices are boundary vertices, then don't
            % do the swap. In the future maybe remove this restriction. The
            % reason would be to avoid very small triangles (poor quality)
            % at the boundary if the boundary is very subdivided.
            continue;
        end

        new_edge_verts_adjtri = MESH_new.verts_tri(new_edge);
        new_edge_verts_valency = cellfun('length', new_edge_verts_adjtri);
        filter_valency = new_edge_verts_valency >= max_valency;
        if any(filter_valency)
            continue;
        end
        
        MESH_trial = MESH_new;
        MESH_trial.tri_timer(longest_edge_adjtri) = 0;

        MESH_trial.tri_verts(longest_edge_adjtri(1),:) = sort([new_edge, longest_edge(1)]);
        MESH_trial.tri_verts(longest_edge_adjtri(2),:) = sort([new_edge, longest_edge(2)]);

        %%%% (REPLACES setdiff operation)
        temp = sort(MESH_trial.verts_tri{longest_edge(1)});
        temp(temp == longest_edge_adjtri(2)) = [];
        MESH_trial.verts_tri{longest_edge(1)} = temp;

        temp = sort(MESH_trial.verts_tri{longest_edge(2)});
        temp(temp == longest_edge_adjtri(1)) = [];
        MESH_trial.verts_tri{longest_edge(2)} = temp;
        %%%%

        %%%% (REPLACES union operation)
        temp = sort([MESH_trial.verts_tri{new_edge(1)}; longest_edge_adjtri]);
        temp_logical = [false; diff(temp) == 0];
        MESH_trial.verts_tri{new_edge(1)} = temp(~temp_logical);

        temp = sort([MESH_trial.verts_tri{new_edge(2)}; longest_edge_adjtri]);
        temp_logical = [false; diff(temp) == 0];
        MESH_trial.verts_tri{new_edge(2)} = temp(~temp_logical);
        %%%%   

        % Check to make sure that new triangles have higher average quality
        % than old ones

        % First check if the quality of the adjacent triangles need to be
        % recomputed
        
        quality_adjtri = zeros(1, 2);
        for jj = 1:2
            if recompute_tri_quality(longest_edge_adjtri(jj))
                verts_ind_temp = MESH_new.tri_verts(longest_edge_adjtri(jj),:);
                verts_temp = MESH_new.verts(:,verts_ind_temp);
    
                edges_temp_length = sqrt(sum([(verts_temp(:,1) - verts_temp(:,2)).^2, ...
                                              (verts_temp(:,1) - verts_temp(:,3)).^2, ...
                                              (verts_temp(:,2) - verts_temp(:,3)).^2], 1));
    
                tri_edge_length(:,longest_edge_adjtri(jj)) = edges_temp_length;
            
                detTxy = verts_temp(1,1)*verts_temp(2,2) - verts_temp(1,2)*verts_temp(2,1) - ...
                         verts_temp(1,1)*verts_temp(2,3) + verts_temp(1,3)*verts_temp(2,1) + ...
                         verts_temp(1,2)*verts_temp(2,3) - verts_temp(1,3)*verts_temp(2,2);
            
                detTyz = verts_temp(2,1)*verts_temp(3,2) - verts_temp(2,2)*verts_temp(3,1) - ...
                         verts_temp(2,1)*verts_temp(3,3) + verts_temp(2,3)*verts_temp(3,1) + ...
                         verts_temp(2,2)*verts_temp(3,3) - verts_temp(2,3)*verts_temp(3,2);
                
                detTzx = verts_temp(1,2)*verts_temp(3,1) - verts_temp(1,1)*verts_temp(3,2) + ...
                         verts_temp(1,1)*verts_temp(3,3) - verts_temp(1,3)*verts_temp(3,1) - ...
                         verts_temp(1,2)*verts_temp(3,3) + verts_temp(1,3)*verts_temp(3,2);
            
                area = 1/2 * sqrt(detTxy^2 + detTyz^2 + detTzx^2);
                perimeter = sum(edges_temp_length);
            
                quality_adjtri(jj) = 36 * area / (sqrt(3) * perimeter^2);
        
                tri_area(longest_edge_adjtri(jj)) = area;
                tri_perim(longest_edge_adjtri(jj)) = perimeter;
                tri_quality(longest_edge_adjtri(jj)) = quality_adjtri(jj);
    
                recompute_tri_quality(longest_edge_adjtri(jj)) = false;
            else
                quality_adjtri(jj) = tri_quality(longest_edge_adjtri(jj));
            end
        end

        % Then compute quality of trial triangles
        quality_adjtri_trial = zeros(1, 2);
        tri_edge_length_trial = zeros(3, 2);
        tri_area_trial = zeros(1, 2);
        tri_perim_trial = zeros(1,2);
        for jj = 1:2
            verts_ind_temp = MESH_trial.tri_verts(longest_edge_adjtri(jj),:);
            verts_temp = MESH_trial.verts(:,verts_ind_temp);

            edges_temp_length = sqrt(sum([(verts_temp(:,1) - verts_temp(:,2)).^2, ...
                                          (verts_temp(:,1) - verts_temp(:,3)).^2, ...
                                          (verts_temp(:,2) - verts_temp(:,3)).^2], 1));
        
            detTxy = verts_temp(1,1)*verts_temp(2,2) - verts_temp(1,2)*verts_temp(2,1) - ...
                     verts_temp(1,1)*verts_temp(2,3) + verts_temp(1,3)*verts_temp(2,1) + ...
                     verts_temp(1,2)*verts_temp(2,3) - verts_temp(1,3)*verts_temp(2,2);
        
            detTyz = verts_temp(2,1)*verts_temp(3,2) - verts_temp(2,2)*verts_temp(3,1) - ...
                     verts_temp(2,1)*verts_temp(3,3) + verts_temp(2,3)*verts_temp(3,1) + ...
                     verts_temp(2,2)*verts_temp(3,3) - verts_temp(2,3)*verts_temp(3,2);
            
            detTzx = verts_temp(1,2)*verts_temp(3,1) - verts_temp(1,1)*verts_temp(3,2) + ...
                     verts_temp(1,1)*verts_temp(3,3) - verts_temp(1,3)*verts_temp(3,1) - ...
                     verts_temp(1,2)*verts_temp(3,3) + verts_temp(1,3)*verts_temp(3,2);
        
            area = 1/2 * sqrt(detTxy^2 + detTyz^2 + detTzx^2);
            perimeter = sum(edges_temp_length);
        
            quality_adjtri_trial(jj) = 36 * area / (sqrt(3) * perimeter^2);
    
            tri_edge_length_trial(:,jj) = edges_temp_length;
            tri_area_trial(jj) = area;
            tri_perim_trial(jj) = perimeter;
        end

        if mean(quality_adjtri_trial) > mean(quality_adjtri)
            MESH_new = MESH_trial;

            tri_quality(longest_edge_adjtri) = quality_adjtri_trial;
            tri_area(longest_edge_adjtri) = tri_area_trial;
            tri_perim(longest_edge_adjtri) = tri_perim_trial;
            tri_edge_length(:,longest_edge_adjtri) = tri_edge_length_trial;

            recompute_tri_quality(longest_edge_adjtri) = false;
        end
        

        % MESH_new = MESH_trial;
        % recompute_tri_quality(longest_edge_adjtri) = true;

        swap_count = swap_count + 1;
    end

end

