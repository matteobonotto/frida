function [G_Br_VI,G_Bz_VI] = fun_G_BrBz_active_VI_stable_fast(n_act, ...
    tri_act, ...
    nodes_act, ...
    ind_act, ...
    keyreg_act, ...
    degree_G_source, ...
    r_point, ...
    z_point)


%%

npt_point = length(r_point);

G_Br_VI = zeros(npt_point,n_act);
G_Bz_VI = zeros(npt_point,n_act);

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
    
    [w_G_act,P_G_act,n_G_act] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source);
        
    % Green matrix
    r_source = P_G_act(:,1);
    z_source = P_G_act(:,2);
    I_source = w_G_act/sum(w_G_act);
    npt_source = length(r_source);
        
    [vec_Br_all,vec_Bz_all] = fun_Green_filament_BrBz_SP_f90(npt_source, ...
        r_source, ...
        z_source, ...
        I_source, ...
        npt_point, ...
        r_point, ...
        z_point, ...
        1,...
        12);
    
    G_Br_VI(:,ii) = vec_Br_all;
    G_Bz_VI(:,ii) = vec_Bz_all;
    

end










































