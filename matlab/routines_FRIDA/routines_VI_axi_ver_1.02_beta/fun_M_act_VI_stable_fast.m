function [M_VI,Area_coil] = fun_M_act_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source_ext, ...
    degree_G_target, ...
    n_G_target, ....
    P_G_target,...
    n_act, ...
    tri_act, ...
    nodes_act, ...
    ind_act, ...
    keyreg_act)


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

[w_G_target_norm,P_G_target_norm,~] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_target);


% Jacobian of linear transformation (area of triangles)
% % edge_1 = P2 - P1;
% % edge_2 = P3 - P2;

edge_1 = nodes(tri(:,2),:) - nodes(tri(:,1),:);
edge_2 = nodes(tri(:,3),:) - nodes(tri(:,2),:);

det_Jac = .5*abs(edge_1(:,1).*edge_2(:,2) - edge_1(:,2).*edge_2(:,1));


%
W_r_target = fun_calc_shape_functions_points(shape_f_norm,P_G_target_norm,N_order);





%%
nn = size(nodes,1);
nt = size(tri,1);

M_VI = zeros(nn,n_act);

Area_coil = zeros(n_act,1);

mu0 = 4*pi*1e-7;

for ii = 1:numel(ind_act)

    tri_act_ii = (tri_act(keyreg_act == ind_act(ii),:));


    %Gauss points for coil ii
    nodes_matrix = [nodes_act(tri_act_ii(:,1),:) ...
        nodes_act(tri_act_ii(:,2),:) ...
        nodes_act(tri_act_ii(:,3),:)];

    P1 = nodes_matrix(:,[1 2]);
    P2 = nodes_matrix(:,[3 4]);
    P3 = nodes_matrix(:,[5 6]);

    [w_G_act,P_G_act,n_G_act] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source_ext);

    Area_coil(ii) = sum(w_G_act);

    % Green matrix
    r_source = P_G_act(:,1);
    z_source = P_G_act(:,2);
    I_source = w_G_act/sum(w_G_act);
    % %     I_source = w_G_act;
    npt_source = length(r_source);

    r_point = P_G_target(:,1);
    z_point = P_G_target(:,2);
    npt_point = length(r_point);

    vec_Aphi_all = fun_Green_filament_Aphi_SP_f90(npt_source, ...
        r_source, ...
        z_source, ...
        I_source, ...
        npt_point, ...
        r_point, ...
        z_point, ...
        1,...
        12);

    % %     figure
    % %     plot(r_point,z_point,'ro'); hold on;
    % %     plot(r_source,z_source,'ok'); hold on;
    % %     title(num2str(ii))

    % Assembly M_loc

    if isfile('fun_assembly_M_VI_mex')
        M_VI = fun_assembly_M_VI_mex(M_VI, ...
            ii, ...
            tri, ...
            n_G_target, ...
            P_G_target, ...
            w_G_target_norm, ...
            vec_Aphi_all, ...
            N_nodes_sf, ...
            W_r_target, ...
            det_Jac);
    else
        M_VI = fun_assembly_M_VI(M_VI, ...
            ii, ...
            tri, ...
            n_G_target, ...
            P_G_target, ...
            w_G_target_norm, ...
            vec_Aphi_all, ...
            N_nodes_sf, ...
            W_r_target, ...
            det_Jac);
    end
    


end

