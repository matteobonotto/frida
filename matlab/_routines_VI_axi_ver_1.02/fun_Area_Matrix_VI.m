function M_A = fun_Area_Matrix_VI(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source)



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

.5*abs(edge_1(1,1).*edge_2(1,2) - edge_1(1,2).*edge_2(1,1)) - det_Jac(1)


%
W_r_source = fun_calc_shape_functions_points(shape_f_norm,P_G_soruce_norm,N_order);


%%
nn = size(nodes,1);

% % M_A = zeros(nn);

nt = size(tri,1);

M_A_index = zeros(nt,N_nodes_sf);
M_A_vals = zeros(nt,N_nodes_sf);

for ii=1:nt
    
    ind_t_pas_ii = ii;
    tri_ii = tri(ind_t_pas_ii,:);
                
    w_G_source_ii = 2*w_G_soruce_norm*det_Jac(ii);
    
    ind_G_ii = (ii-1)*n_G_source+1:ii*n_G_source;
    P_G_source_ii = P_G_source(ind_G_ii,:);
    r_G_source_ii = P_G_source_ii(:,1);

% %     if ismember(5974,tri_ii)
% %     aa = 0;
% %     end
% %     
% %     P1 = nodes(tri_ii(1),:);
% %     P2 = nodes(tri_ii(2),:);
% %     P3 = nodes(tri_ii(3),:);
% % 
% %     [w_G_soruce,P_G_soruce,~] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source);
% %     [shape_f_ii] = fun_shape_functions_stable(tri_ii,nodes,N_order);
% %     W_r_source_ii = fun_calc_shape_functions_points(shape_f_ii,P_G_soruce,N_order);
        
    M_A_loc = zeros(N_nodes_sf,1);
    for hh = 1:N_nodes_sf
        
% %         figure
% %         plot3(P_G_soruce_norm(:,1),P_G_soruce_norm(:,2),W_r_source(:,hh),'o')
% %         hold on;
% %         plot(nodes_norm(:,1),nodes_norm(:,2), 'o')
% %         
% %         k_r = [2 0 1 1 0 0];
% %         k_z = [0 2 1 0 1 0];
% %         cc_shp_ii = reshape(shape_f_norm,6,6).';
% %         
% %         res = 0;
% %         for kk = 1:6
% %             res = res + 2*det_Jac(ii)*cc_shp_ii(hh,kk)*factorial(k_r(kk))*factorial(k_z(kk))/factorial(2+k_r(kk)+k_z(kk));
% %         end
% %         
% %         W_r_source(:,hh)'*w_G_source_ii
% %         W_r_source(:,hh)'*w_G_soruce
% %         W_r_source_ii(:,hh)'*w_G_soruce
% %         
% %         xx = linspace(0,1,100)
% %         figure
% %         plot(xx,2*xx.^2-3*xx+1,'.')
% %         fun = @(x) x.^2 - 3*x +1;
% %         q = integral(fun,0,1)
% %         
% %         fun = @(x,y) cc_shp_ii(hh,1)*x.^2 + ...
% %             cc_shp_ii(hh,2)*y.^2 + ...
% %             cc_shp_ii(hh,3)*x.*y + ...
% %             cc_shp_ii(hh,4)*x + ...
% %             cc_shp_ii(hh,5)*y + ...
% %             cc_shp_ii(hh,6)
% %         xmin = 0;
% %         xmax = 1;
% %         ymin = 0;
% %         ymax = @(x) 1-x;
% %         Int_matlab = integral2(fun,xmin,xmax,ymin,ymax,'AbsTol', 0,'RelTol',eps);

        
        
        M_A_loc(hh) = (2*pi*r_G_source_ii.*W_r_source(:,hh))'*w_G_source_ii;
    end
    
    M_A_index(ii,:)  = tri_ii;
    M_A_vals(ii,:) = M_A_loc.';
      
% %     M_A_loc = zeros(3*N_order);
% %     for hh = 1:N_nodes_sf
% %         for kk = 1:N_nodes_sf
% % 
% %             temp_int = W_r_source(:,hh).*W_r_source(:,kk);
% %                         
% %             M_A_loc(hh,kk) = temp_int'*w_G_source_ii;
% %             
% %         end
% %     end
% %     
% %     M_A(tri_ii,tri_ii) = M_A(tri_ii,tri_ii) + M_A_loc;
    
end


M_A = sparse(M_A_index(:),M_A_index(:),M_A_vals(:));




