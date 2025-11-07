function G_VI = fun_assemby_G_VI(tri, ...
    nn, ...
    n_G_source, ...
    Green_Mat_Gauss, ...
    det_Jac, ...
    w_G_soruce_norm, ...
    W_r_source)

%%

nt = size(tri,1);
npt_point = size(Green_Mat_Gauss,1);

N_nodes_sf = size(W_r_source,2);

G_VI = zeros(npt_point,nn);

for ii=1:nt
    
    ind_t_pas_ii = ii;
    tri_ii = tri(ind_t_pas_ii,:);
    
    ind_G_ii = (ii-1)*n_G_source+1:ii*n_G_source;
    
    w_G_source_ii = 2*w_G_soruce_norm*det_Jac(ii);
        
    flux_ii = Green_Mat_Gauss(:,ind_G_ii);

    
    %
    G_loc = zeros(npt_point,N_nodes_sf);
    for hh = 1:N_nodes_sf
        
        temp_product = (flux_ii.*repmat(W_r_source(:,hh).',npt_point,1));
        G_loc(:,hh) = temp_product*(w_G_source_ii);
        
    end
    
    G_VI(:,tri_ii) = G_VI(:,tri_ii) + G_loc;
    
    
end


