function L_VI = fun_assemby_L_VI(tri, ...
    nn, ...
    n_G_source, ...
    n_G_target, ....
    P_G_target, ...
    Green_Mat_Gauss_Aphi, ...
    det_Jac, ...
    w_G_soruce_norm, ...
    W_r_source, ...
    w_G_target_norm, ...
    W_r_target)
    

%%

nt = size(tri,1);

N_nodes_sf = size(W_r_source,2);

L_VI = zeros(nn);

for ii=1:nt
    
    ind_t_pas_ii = ii;
    tri_ii = tri(ind_t_pas_ii,:);
        
    ind_G_ii = (ii-1)*n_G_source+1:ii*n_G_source;
        
    w_G_source_ii = 2*w_G_soruce_norm*det_Jac(ii);
% %     w_G_source_ii = 2*w_G_soruce_norm;
    
% %     J_source_ii = 1/sum(w_G_source_ii);
    
% %     P_G_source_ii = P_G_source(ind_G_ii,:);
    
    for jj=1:nt
    
        ind_t_pas_jj = jj;
        tri_jj = tri(ind_t_pas_jj,:);
        
        ind_G_jj = (jj-1)*n_G_target+1:jj*n_G_target;
                
        w_G_target_jj = 2*w_G_target_norm*det_Jac(jj);
% %         w_G_target_jj = 2*w_G_target_norm;

% %         J_target_jj = 1/sum(w_G_target_jj);

        
        P_G_target_jj = P_G_target(ind_G_jj,:);
        r_G_target_jj = P_G_target_jj(:,1);
        
        L_loc = zeros(N_nodes_sf);
        
        Aphi_ii_jj = Green_Mat_Gauss_Aphi(ind_G_jj,ind_G_ii);
        
        
        for hh = 1:N_nodes_sf

            temp_product = (2*pi*Aphi_ii_jj.*repmat(W_r_source(:,hh).',n_G_target,1));
            temp_int_hh = temp_product*(w_G_source_ii);           
            
            for kk = 1:N_nodes_sf

                temp_int = (r_G_target_jj.*temp_int_hh.*W_r_target(:,kk)).'*w_G_target_jj;

                L_loc(hh,kk) = temp_int;
                
            end
        end
                
        L_VI(tri_ii,tri_jj) = L_VI(tri_ii,tri_jj) + L_loc;

    
    end
    
end








