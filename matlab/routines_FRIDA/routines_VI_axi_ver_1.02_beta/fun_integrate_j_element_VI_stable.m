function I_element = fun_integrate_j_element_VI_stable(tri, ...
    nodes, ...
    j_VI, ...
    N_order, ...
    degree_G_source, ...
    ind_element, ...
    type)



%% Quantities on normalized triangle

[shape_f_norm,nodes_norm,N_nodes_sf] = fun_shape_functions_norm(N_order); %#

if N_order == 2
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
% % W_r_source = fun_calc_shape_functions_points(shape_f_norm,P_G_soruce_norm,N_order);


%%

I_element = zeros(size(ind_element));

for jj = 1:length(ind_element)
    
    ind_t_sel = find(type == ind_element(jj));
% %     tri_sel = tri(ind_t_sel,:);
    
    I_element_jj = 0;
    
    nt_sel = length(ind_t_sel);
    
    for ii = 1:nt_sel
        
        ind_t_pas_ii = ind_t_sel(ii);
        tri_ii = tri(ind_t_pas_ii,:);
        
        w_G_source_ii = 2*w_G_soruce_norm*det_Jac(ind_t_pas_ii);
                
        j_VI_ii = j_VI(tri_ii);

        j_VI_G_ii = fun_calc_function_on_points(shape_f_norm,P_G_soruce_norm,j_VI_ii,N_order);
        
% %         figure
% %         plot3(nodes_norm(:,1),nodes_norm(:,2),j_VI_ii,'o')
% %         hold on
% %         plot3(P_G_soruce_norm(:,1),P_G_soruce_norm(:,2),j_VI_G_ii,'o')
        
        I_element_ii = w_G_source_ii'*j_VI_G_ii;
        
        %
        I_element_jj = I_element_jj + I_element_ii;
        
% %         pause
        
    end
    
    I_element(jj) = I_element_jj;
    
end


% % shape_f_vec = shape_f_norm;
% % P_target = P_G_soruce_norm;
% % npt_shape_functions = sqrt(numel(shape_f_vec));
% % n_target = size(P_target,1);
% % 
% % cc_shp_ii_norm = reshape(shape_f_vec,npt_shape_functions,npt_shape_functions).';
% % 
% % fac_geo_norm = [P_target.^2 ...
% %     P_target(:,1).*P_target(:,2) ...
% %     P_target ...
% %     ones(n_target,1)].';
% % 
% % w_target_norm = (cc_shp_ii_norm*fac_geo_norm);
% % 
% % cc_shp_ii_norm(4,:)*fac_geo_norm(:,1)
% % cc_shp_ii_norm(5,:)*fac_geo_norm(:,1)
% % 
% % 
% % shape_f_vec = shape_f_ii;
% % P_target = P_G_source_ii;
% % npt_shape_functions = sqrt(numel(shape_f_vec));
% % n_target = size(P_target,1);
% % 
% % cc_shp_ii = reshape(shape_f_vec,npt_shape_functions,npt_shape_functions).';
% % 
% % fac_geo = [P_target.^2 ...
% %     P_target(:,1).*P_target(:,2) ...
% %     P_target ...
% %     ones(n_target,1)].';
% % 
% % w_target = (cc_shp_ii*fac_geo);
% % cc_shp_ii(4,:)*fac_geo(:,1)
% % cc_shp_ii(5,:)*fac_geo(:,1)
% % 
% % 
% % nodes_1 = nodes_norm;
% % [aa,bb] = fun_shape_functions([1 2 3 4 5 6],nodes_1,N_order)

% % figure
% % triplot(tri(1,1:3),nodes(:,1),nodes(:,2), 'LineWidth',2)
% % hold on
% % plot(nodes(tri(1,1),1),nodes(tri(1,1),2),'o', 'LineWidth',2)
% % plot(nodes(tri(1,2),1),nodes(tri(1,2),2),'o', 'LineWidth',2)
% % plot(nodes(tri(1,3),1),nodes(tri(1,3),2),'o', 'LineWidth',2)
% % plot(nodes(tri(1,4),1),nodes(tri(1,4),2),'o', 'LineWidth',2)
% % plot(nodes(tri(1,5),1),nodes(tri(1,5),2),'o', 'LineWidth',2)
% % 
% % 
% % for ii = 1:6
% %     
% %     figure
% %     plot(nodes_norm([1:3 1],1),nodes_norm([1:3 1],2),'k'); hold on
% %     plot(nodes_norm(ii,1),nodes_norm(ii,2),'o', 'LineWidth',2)
% %     
% %     plot3(P_G_soruce_norm(:,1),P_G_soruce_norm(:,2),W_r_source(:,ii),'o')
% %     view(3)
% %     pause
% %     
% % end
% % 
% % for ii = 1:6
% %     
% %     figure
% %     triplot(tri(1,1:3),nodes(:,1),nodes(:,2), 'LineWidth',2)
% %     hold on
% %     plot(nodes(tri(1,ii),1),nodes(tri(1,ii),2),'o', 'LineWidth',2)
% %     
% %     plot3(P_G_source_ii(:,1),P_G_source_ii(:,2),W_r_hh(:,ii),'o')
% %     view(3)
% %     pause
% %     
% % end











