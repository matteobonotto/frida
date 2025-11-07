function [Psi_Lim,RZ_Lim] = fun_find_limiter_point(Axis,curve,psi_curve,grad_curve)


cenro_curve = sum(curve)/size(curve,1);
curve_0 = curve - sum(curve)/size(curve,1);
theta = atan2(curve_0(:,2),curve_0(:,1));
rho = sqrt(curve_0(:,1).^2 + curve_0(:,2).^2);

[theta_sort,ind] = sort(theta);
psi_curve_sort = psi_curve(ind);
rho_curve_sort = rho(ind);

theta_sort = [theta_sort(end)-2*pi; theta_sort; theta_sort(1)+2*pi;];
psi_curve_sort = psi_curve_sort([end 1:end 1]);
rho_curve_sort = rho_curve_sort([end 1:end 1]);

%%

prod_scal = sum((curve - Axis).*grad_curve,2);


% % figure
% % plot3(curve(:,1),curve(:,2),psi_curve,'o'); hold on;
% % plot3(curve(prod_scal>0,1),curve(prod_scal>0,2),psi_curve(prod_scal>0),'*'); hold on;

% % plot(Axis(1),Axis(2),'o'); hold on;
% % plot(curve(prod_scal>0,1),curve(prod_scal>0,2),'*'); hold on;
% % quiver(curve(:,1),curve(:,2),grad_curve(:,1),grad_curve(:,2),'k'); hold on;

prod_scal_sort = prod_scal(ind);
prod_scal_sort = prod_scal_sort([end 1:end 1]);

ind_not = find(prod_scal_sort > 0);
theta_not_min = min(theta_sort(ind_not));
theta_not_max = max(theta_sort(ind_not));

%%
n_piece = 10;
theta_new = linspace(-pi,pi,3*n_piece+1).';
psi_curve_new = interp1(theta_sort,psi_curve_sort,theta_new);

theta_new_last = theta_new(end);
psi_curve_new_last = psi_curve_new(end);

theta_new = theta_new(1:end-1);
psi_curve_new = psi_curve_new(1:end-1);

theta_rshp = reshape(theta_new,3,numel(theta_new)/3).';
Psi_rshp = reshape(psi_curve_new,3,numel(theta_new)/3).';

theta_rshp = [theta_rshp theta_rshp([2:end 1],1)].';
Psi_rshp = [Psi_rshp Psi_rshp([2:end 1],1)].';

theta_rshp(end,end) = theta_new_last;
Psi_rshp(end,end) = psi_curve_new_last;


%%
% % figure
% % plot(theta_new,psi_curve_new)
% % hold on

SP_candidates_theta = NaN(n_piece,2);
SP_candidates_psi = NaN(n_piece,2);

for ii = 1:n_piece
    coeffs = fun_spline_1D_coeffs_f_vals(Psi_rshp(:,ii),theta_rshp(:,ii));
% %     theta_ii = linspace(theta_rshp(1,ii),theta_rshp(end,ii),20).';
% %     ft = fun_eval_spline_1D(coeffs,theta_ii);
% %     pause
% %     plot(theta_ii,ft,'o')
    
    a = 3*coeffs(4);
    b = 2*coeffs(3);
    c = coeffs(2);
    delta = b^2 - 4*a*c;
    x_12 = [-b+sqrt(delta) -b-sqrt(delta)]/2/a;
    
    x_12(find(imag(x_12))) = NaN;
    
    if isempty(theta_not_min)
        is_in = x_12(:) >= theta_rshp(1,ii) & x_12(:) < theta_rshp(end,ii);
    else
        is_in = x_12(:) >= theta_rshp(1,ii) & x_12(:) < theta_rshp(end,ii) ...
            & ~(x_12(:) > theta_not_min & x_12(:) < theta_not_max);
    end
    
    if any(is_in)
        SP_candidates_theta(ii,1:numel(find(is_in))) = x_12(is_in);
        % %         plot(x_12(is_in),fun_eval_spline_1D(coeffs,x_12(is_in)),'*r', 'LineWidth',2)
        SP_candidates_psi(ii,1:numel(find(is_in))) = fun_eval_spline_1D(coeffs,x_12(is_in));
    end
    % %     pause
    
end

theta_candidates = SP_candidates_theta(:);
theta_candidates = theta_candidates(~isnan(theta_candidates));

psi_candidates = SP_candidates_psi(:);
psi_candidates = psi_candidates(~isnan(psi_candidates));


[Psi_Lim,ind_Lim] = max(psi_candidates);
rho_lim = interp1(theta_sort,rho_curve_sort,theta_candidates(ind_Lim));
RZ_Lim = rho_lim*[cos(theta_candidates(ind_Lim)) sin(theta_candidates(ind_Lim))] + sum(curve)/size(curve,1);


% % figure
% % plot3(curve(:,1),curve(:,2),psi_curve,'o'); hold on;
% % plot3(curve(:,1),curve(:,2),psi_curve,'o'); hold on;
% % plot3(RZ_Lim(1),RZ_Lim(2),Psi_Lim,'*', 'LineWidth',2);
% % 
% % 
% % plot(RZ_Lim(1),RZ_Lim(2),'*', 'LineWidth',2);























