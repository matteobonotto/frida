function M_A = fun_vec_Area_element_VI_stable(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    ind_element, ...
    type) %#codegen


%% Quantities on normalized triangle

[shape_f_norm,nodes_norm,N_nodes_sf] = fun_shape_functions_norm(N_order); %#

if N_order > 1
    midpoints = .5*(nodes(tri(1,1:3),:) + nodes(tri(1,[2 3 1]),:));
    
    ind_map = zeros(N_nodes_sf,1);
    ind_map(1:3) = 1:3;
    if norm(midpoints - nodes(tri(1,4:end),:)) > 1e-5
        
        vec = nodes(tri(1,4),:);
        temp_mat = midpoints - repmat(vec,size(midpoints,1),1);
        
        [~,ind_min] = min(sqrt(temp_mat(:,1).^2 + temp_mat(:,2).^2));
        
        ind_map(4) = ind_min+3;
        ind_map(5:end) = ind_map(4)+1:ind_map(4)+(N_nodes_sf-3-1);
        ind_map(ind_map>N_nodes_sf) = ind_map(ind_map>N_nodes_sf) - 3;
        
        temp_shape_f_norm = shape_f_norm(:);
        temp_shape_f_norm = reshape(temp_shape_f_norm,N_nodes_sf,N_nodes_sf).';
        
        temp_shape_f_norm = temp_shape_f_norm(ind_map,:);
        shape_f_norm = reshape(temp_shape_f_norm.',1,N_nodes_sf^2);
        
        nodes_norm = nodes_norm(ind_map,:);
        
    end
end

P1 = nodes_norm(1,:);
P2 = nodes_norm(2,:);
P3 = nodes_norm(3,:);

[w_G_soruce_norm,P_G_soruce_norm,~] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source);


% Jacobian of linear transformation (area of triangles)
% % edge_1 = P2 - P1;
% % edge_2 = P3 - P2;

edge_1 = nodes(tri(:,2),:) - nodes(tri(:,1),:);
edge_2 = nodes(tri(:,3),:) - nodes(tri(:,2),:);

det_Jac = .5*abs(edge_1(:,1).*edge_2(:,2) - edge_1(:,2).*edge_2(:,1));


%
W_r_source = fun_calc_shape_functions_points(shape_f_norm,P_G_soruce_norm,N_order);


%%
nn = size(nodes,1);
M_A = zeros(nn,length(ind_element));

for jj = 1:length(ind_element)
    
    ind_t_sel = find(type == ind_element(jj));
    
    ind_n_sel = unique(tri(ind_t_sel,:));
    
    M_A_jj = zeros(nn,1);
    nn_sel =length(ind_n_sel);
    
    for ii=1:nn_sel
        
        ind_t_base_source_ii = MatInd_nodes_tri(ind_n_sel(ii),:);
        ind_t_base_source_ii = ind_t_base_source_ii(ind_t_base_source_ii ~= 0);
        
        nt_base_source = length(ind_t_base_source_ii);
        
        ind_n_base_source_ii_Mat = tri(ind_t_base_source_ii,:);
        
        ind_G_base_source_ii = zeros(length(ind_t_base_source_ii),n_G_source);
        for hh = 1:length(ind_t_base_source_ii)
            ind_G_base_source_ii(hh,:) = (ind_t_base_source_ii(hh)-1)*n_G_source+1:ind_t_base_source_ii(hh)*n_G_source;
        end
        ind_G_base_source_ii = ind_G_base_source_ii.';
        
        for hh = 1:nt_base_source
            
            ind_sel_hh = (hh-1)*n_G_source+1:hh*n_G_source;
            ind_nodes_sel_hh = (ind_n_base_source_ii_Mat(hh,:) == ind_n_sel(ii));
            
% %             if sum(ind_nodes_sel_hh) == 0
% %                 res_hh = 0;
% %             else
                w_G_base_source_ii_hh = 2*w_G_soruce_norm.*det_Jac(ind_t_base_source_ii(hh));
                res_hh = W_r_source(:,ind_nodes_sel_hh).'*w_G_base_source_ii_hh;
% %             end
            
            M_A_jj(ind_n_sel(ii)) = M_A_jj(ind_n_sel(ii)) + res_hh;
            
        end
        
    end
    
    M_A(:,jj) = M_A_jj;
    
    
end

end














