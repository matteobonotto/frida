function [ind_t_on_curve] = fun_ind_t_on_separatrix(tri,ind_n_vertices_plasma,ind_n_vertices_not_plasma)

nt = size(tri,1);
ind_t_on_curve = zeros(nt,1);
tri_ord_1 = tri(:,1:3);

for jj = 1:nt % parfor nella versione mexata
    
    tri_jj = tri_ord_1(jj,:);
    
    test_in = ismember(tri_jj',ind_n_vertices_plasma);
    test_out = ismember(tri_jj',ind_n_vertices_not_plasma);
    
    if sum(test_in)>0 && sum(test_out)>0
        ind_t_on_curve(jj) = 1;
    end
    
end

ind_t_on_curve = find(ind_t_on_curve);