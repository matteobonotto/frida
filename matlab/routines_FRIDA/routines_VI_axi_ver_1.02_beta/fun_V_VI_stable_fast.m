function V_VI = fun_V_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source)



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
W_r_source = fun_calc_shape_functions_points(shape_f_norm,P_G_soruce_norm,N_order);


%%
nn = size(nodes,1);

V_VI = zeros(nn);

nt = size(tri,1);

for ii=1:nt
    
    ind_t_pas_ii = ii;
    tri_ii = tri(ind_t_pas_ii,:);
        
    ind_G_ii = (ii-1)*n_G_source+1:ii*n_G_source;
    
    ind_G_ii = ind_G_ii(:);
    
    w_G_source_ii = 2*w_G_soruce_norm*det_Jac(ii);

    % %     nodes_ii = nodes(tri_ii,:);
    % %     r1 = nodes_ii(1,1);
    % %     r21 = nodes_ii(2,1) - nodes_ii(1,1);
    % %     r31 = nodes_ii(3,1) - nodes_ii(1,1);
    
    
    P_G_source_ii = P_G_source(ind_G_ii,:);
    r_G_source_ii = P_G_source_ii(:,1);
    
    U_loc = zeros(3*N_order);
    
    for hh = 1:N_nodes_sf
        for kk = 1:N_nodes_sf
            
            temp_int = W_r_source(:,hh).*W_r_source(:,kk);
% %             temp_int = W_r_source(:,hh).*W_r_source(:,kk);
            
            U_loc(hh,kk) = temp_int'*w_G_source_ii;
            
            % Analytical computation
            % %             m = [2 0 1 1 0 0];
            % %             n = [0 2 1 0 1 0];
            % %             cc_shp_ii = reshape(shape_f_norm,6,6).';
            % %
            % %             aa_1 = cc_shp_ii(hh,:);
            % %             aa_2 = cc_shp_ii(kk,:);
            % %             res_int = 0;
            % %             for ss = 1:6
            % %                 for tt = 1:6
            % %                     exp_r = m(ss)+m(tt);
            % %                     exp_z = n(ss)+n(tt);
            % %
            % %                     I_1 = r1*fun_integrate_poly_on_tri(exp_r,exp_z);
            % %                     I_2 = r21*fun_integrate_poly_on_tri(exp_r+1,exp_z);
            % %                     I_3 = r31*fun_integrate_poly_on_tri(exp_r,exp_z+2);
            % %
            % %                     temp_res_int = 2*pi*aa_1(ss)*aa_2(tt)*2*det_Jac(ii)*(I_1+I_2+I_3);
            % %
            % %
            % %                     res_int = res_int + temp_res_int;
            % %
            % %                 end
            % %             end
            
        end
    end
    
    V_VI(tri_ii,tri_ii) = V_VI(tri_ii,tri_ii) + U_loc;
    
end













