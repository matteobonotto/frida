function [ind_t_B] = fun_indTriOnEdge(tri,ind_n_B)

nt = size(tri,1);

tmp_ind_t_bound = zeros(nt,1);

for ii=1:nt % parfor in the mex version

    tri_ii = tri(ii,:)';
    
    tmp = ismember(tri_ii,ind_n_B);
    if sum(tmp) > 1
        tmp_ind_t_bound(ii) = 1;
    end

end

%
ind_t_B = find(tmp_ind_t_bound);