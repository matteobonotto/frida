function D_VI = fun_D_VI(tri,type_pas,ind_t_map)

nn = max(tri(:));

D_VI = zeros(nn,size(ind_t_map,1));

for ii = 1:size(ind_t_map,1)
    
    if ind_t_map(ii,1) == ind_t_map(ii,2)
        ind_sel = find((type_pas == (ind_t_map(ii,1):ind_t_map(ii,2))));
    else
        ind_sel = find(sum(type_pas.' == (ind_t_map(ii,1):ind_t_map(ii,2)).'))';
    end

    D_VI(unique(tri(ind_sel,:)),ii) = 1;
    
end