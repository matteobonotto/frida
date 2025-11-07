function w_target = fun_calc_shape_functions_points(shape_f_vec,P_target,N_order)

npt_shape_functions = sqrt(numel(shape_f_vec));

cc_shp_ii = reshape(shape_f_vec,npt_shape_functions,npt_shape_functions).';

n_target = size(P_target,1);

if N_order == 1
    fac_geo = [P_target ones(n_target,1)].';
    
elseif N_order == 2
    fac_geo = [P_target.^2 ...
        P_target(:,1).*P_target(:,2) ...
        P_target ...
        ones(n_target,1)].';
    
elseif N_order == 3
    fac_geo = [P_target.^3 ...
        P_target(:,1).^2.*P_target(:,2) ...
        P_target(:,1).*P_target(:,2).^2 ...
        P_target.^2 ...
        P_target(:,1).*P_target(:,2) ...
        P_target ...
        ones(n_target,1)].';
    
else
    
    error ('N_order too high, max N_order= 3')
    
end

w_target = (cc_shp_ii*fac_geo).'; 

% % figure
% % plot3(P_target(:,1),P_target(:,2),w_target,'o')

























