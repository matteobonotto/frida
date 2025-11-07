function R_VI = fun_R_VI(tri, ...
    nn,...
    N_order, ...
    n_G_source, ...
    P_G_source,...
    w_G_source,...
    shape_f, ...
    eta)

%%

R_VI = zeros(nn);

nt = size(tri,1);

for ii=1:nt
    
    ind_t_pas_ii = ii;
    tri_ii = tri(ind_t_pas_ii,:);
        
    ind_G_ii = (ii-1)*n_G_source+1:ii*n_G_source;
    
    ind_G_ii = ind_G_ii(:);
    
    w_G_source_ii = w_G_source(ind_G_ii);
    P_G_source_ii = P_G_source(ind_G_ii,:);
    r_G_source_ii = P_G_source_ii(:,1);
    
    shape_f_ii = shape_f(ind_t_pas_ii,:);
    
    W_r_hh = fun_calc_shape_functions_points(shape_f_ii,P_G_source_ii,N_order);
    
    eta_ii = eta(ii);
    
    R_loc = zeros(3*N_order);
    
    for hh = 1:3*N_order
        for kk = 1:3*N_order

            temp_int = 2*pi*r_G_source_ii*eta_ii.*W_r_hh(:,hh).*W_r_hh(:,kk);
            
            R_loc(hh,kk) = temp_int'*w_G_source_ii;
            
        end
    end
    
    R_VI(tri_ii,tri_ii) = R_VI(tri_ii,tri_ii) + R_loc;
    
end