function M_VI = fun_assembly_M_VI(M_VI, ...
    ii, ...
    tri, ...
    n_G_target, ...
    P_G_target, ...
    w_G_target_norm, ...
    vec_Aphi_all, ...
    N_nodes_sf, ...
    W_r_target, ...
    det_Jac)

%%

nt = size(tri,1);
for jj=1:nt
    
    ind_t_pas_jj = jj;
    tri_jj = tri(ind_t_pas_jj,:);
    
    ind_G_jj = (jj-1)*n_G_target+1:jj*n_G_target;
    
    w_G_target_jj = 2*w_G_target_norm*det_Jac(jj);
% %     w_G_target_jj = 2*w_G_target_norm;

% %     J_source_jj = 1/sum(w_G_target_jj);
    
    P_G_target_jj = P_G_target(ind_G_jj,:);
    r_G_target_jj = P_G_target_jj(:,1);
    
    
    Aphi_ii_jj = vec_Aphi_all(ind_G_jj);
    
    M_loc = zeros(N_nodes_sf,1);
    for kk = 1:N_nodes_sf
        
        temp_int = (2*pi*r_G_target_jj.*Aphi_ii_jj.*W_r_target(:,kk)).'*w_G_target_jj;
        
        M_loc(kk) = temp_int;
        
    end
    
    M_VI(tri_jj,ii) = M_VI(tri_jj,ii) + M_loc;
    
end