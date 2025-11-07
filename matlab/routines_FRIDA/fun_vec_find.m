function res = fun_vec_find(ind_1,ind_2)

if size(ind_1,1) < size(ind_1,2)
    ind_1 = ind_1.';
end

if length(ind_2) == 1
    res = find(ind_1 == ind_2);
else
    if size(ind_2,1) < size(ind_2,2)
        ind_2 = ind_2.';
    end
    
    res = find(sum(ind_1.' == ind_2)).';
    
end



