function L_VI =fun_L_VI_stable_symmetric_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    degree_G_target, ...
    n_G_target, ....
    P_G_target)


%%
nt = size(tri,1);
nn = size(nodes,1);

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
[w_G_target_norm,P_G_target_norm,~] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_target);


% Jacobian of linear transformation (area of triangles)
% % edge_1 = P2 - P1;
% % edge_2 = P3 - P2;

edge_1 = nodes(tri(:,2),:) - nodes(tri(:,1),:);
edge_2 = nodes(tri(:,3),:) - nodes(tri(:,2),:);

det_Jac = .5*abs(edge_1(:,1).*edge_2(:,2) - edge_1(:,2).*edge_2(:,1));


%
W_r_source = fun_calc_shape_functions_points(shape_f_norm,P_G_soruce_norm,N_order);
W_r_target = fun_calc_shape_functions_points(shape_f_norm,P_G_target_norm,N_order);



%% Self inductance terms (triangles on themselves)

mu0 = 4*pi*1e-7;


Green_Mat_Gauss_Aphi = zeros(n_G_target,n_G_source*nt);
% % tic
for ii = 1:nt

    ind_G_source = (ii-1)*n_G_source+1:ii*n_G_source;
    ind_G_target = (ii-1)*n_G_target+1:ii*n_G_target;

    for jj=1:length(ind_G_source)
        % source (one Gauss point)
        r_source = P_G_source(ind_G_source(jj),1);
        z_source = P_G_source(ind_G_source(jj),2);
        I_source = 1;
        npt_source = 1;


        % source (other Gauss points of the same triangles)
        r_point = P_G_target(ind_G_target,1);
        z_point = P_G_target(ind_G_target,2);
        npt_point = length(r_point);


        vec_Aphi_all = fun_Green_filament_Aphi_SP_f90(npt_source, ...
            r_source, ...
            z_source, ...
            I_source, ...
            npt_point, ...
            r_point, ...
            z_point, ...
            0,...
            12);

        Green_Mat_Gauss_Aphi(:,ind_G_source(jj)) = vec_Aphi_all;

    end

end
% % toc



% % tic
L_VI_symm = fun_assemby_L_self_VI(tri, ...
    nn, ...
    n_G_source, ...
    n_G_target, ....
    P_G_target, ...
    Green_Mat_Gauss_Aphi, ...
    det_Jac, ...
    w_G_soruce_norm, ...
    W_r_source, ...
    w_G_target_norm, ...
    W_r_target);
% % toc

clear Green_Mat_Gauss_Aphi


%% Mutual inductance terms

r_point = P_G_target(:,1);
z_point = P_G_target(:,2);
npt_point = length(r_point);
Green_Mat_Gauss_Aphi = zeros(n_G_target*nt,n_G_target*nt);

% % tic
for ii = 1:nt

    ind_G = (ii-1)*n_G_target+1:ii*n_G_target;

    for jj=1:length(ind_G)
        % source (one Gauss point)
        r_source = P_G_target(ind_G(jj),1);
        z_source = P_G_target(ind_G(jj),2);
        I_source = 1;
        npt_source = 1;

        vec_Aphi_all = fun_Green_filament_Aphi_SP_f90(npt_source, ...
            r_source, ...
            z_source, ...
            I_source, ...
            npt_point, ...
            r_point, ...
            z_point, ...
            1,...
            12);

        vec_Aphi_all(ind_G) = 0;

        Green_Mat_Gauss_Aphi(:,ind_G(jj)) = vec_Aphi_all;

    end

end
% % toc


if isfile('fun_assemby_L_VI_mex')
    L_VI_mutual = fun_assemby_L_VI_mex(tri, ...
        nn, ...
        n_G_target, ...
        n_G_target, ....
        P_G_target, ...
        Green_Mat_Gauss_Aphi, ...
        det_Jac, ...
        w_G_target_norm, ...
        W_r_target, ...
        w_G_target_norm, ...
        W_r_target);
else
    L_VI_mutual = fun_assemby_L_VI(tri, ...
        nn, ...
        n_G_target, ...
        n_G_target, ....
        P_G_target, ...
        Green_Mat_Gauss_Aphi, ...
        det_Jac, ...
        w_G_target_norm, ...
        W_r_target, ...
        w_G_target_norm, ...
        W_r_target);
end

clear Green_Mat_Gauss_Aphi


%%



L_VI = L_VI_symm + L_VI_mutual;



%%








































