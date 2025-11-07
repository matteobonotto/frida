function [Psi_Lim,RZ_Lim] = fun_find_limiter_point_v3(curve,vers_n,psi_curve,grad_curve)


vers_n = vers_n./fun_vecnorm(vers_n);

[curve_sort,ind] = fun_ordinapunti(curve);
psi_curve_sort = psi_curve(ind);

prod_scal = sum(vers_n.*grad_curve,2);
prod_scal_sort = prod_scal(ind,:);

if norm(curve_sort(1,:) - curve_sort(end,:)) > 1e-6
    curve_sort = [curve_sort; curve_sort(1,:)];
    psi_curve_sort = [psi_curve_sort; psi_curve_sort(1)];
    prod_scal_sort = [prod_scal_sort; prod_scal_sort(1)];
end

curve_sort = [curve_sort; curve_sort(2,:)];
psi_curve_sort = [psi_curve_sort; psi_curve_sort(1)];
prod_scal_sort = [prod_scal_sort; prod_scal_sort(1)];

s_coo = [0; cumsum(fun_vecnorm(curve_sort(2:end,:)-curve_sort(1:end-1,:)))];



%% strategy #2

[s_star,fs_star] = fun_find_SP_1D(s_coo,psi_curve_sort);

if all(prod_scal_sort<0)
    [Psi_Lim,ind_Psi_Lim] = max(fs_star);
    s_star_Lim = s_star(ind_Psi_Lim);
    
else
    prod_scal_sort_s_star = interp1(s_coo,prod_scal_sort,s_star);
    ind_fs_star_pos_prodscal = find(prod_scal_sort_s_star<0);
    [Psi_Lim,ind_Psi_Lim] = max(fs_star(ind_fs_star_pos_prodscal));
    s_star_Lim = s_star(ind_fs_star_pos_prodscal(ind_Psi_Lim));
    
end

RZ_Lim = interp1(s_coo,curve_sort,s_star_Lim);

% % figure
% % plot(s_coo,psi_curve_sort,'ko')
% % hold on;
% % plot(s_star,fs_star,'*r')
% % plot(s_coo(prod_scal_sort>=0),psi_curve_sort(prod_scal_sort>=0),'c*')
% % plot(s_star_Lim,Psi_Lim,'dg')
% % 
% % figure
% % plot3(curve(:,1),curve(:,2),psi_curve,'o')
% % figure
% % plot(curve(:,1),psi_curve,'o')
% % 
% % aa = 0;

end


function OUT = vec_linspace(vec_start,vec_end,nstep)

OUT = zeros(size(vec_end,1),nstep);
for ii=1:size(vec_end,1)
    OUT(ii,:) = linspace(vec_start(ii), vec_end(ii), nstep);
end

end














