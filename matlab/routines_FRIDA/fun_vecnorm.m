function [vec_norm]= fun_vecnorm(vector)

if size(vector,1) == 1 && size(vector,2) == 2
    
    vec_norm = norm(vector);

else
    
    if size(vector,1) < size(vector,2)
        vector = vector';
    end
    
    vec_norm = sqrt(vector(:,1).^2 + vector(:,2).^2);
    
end