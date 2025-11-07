function [Psi_Lim,RZ_Lim] = fun_find_limiter_point_v2(curve,vers_n,psi_curve,grad_curve)


curve_0 = curve - sum(curve)/size(curve,1);
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

% % % % figure
% % % % plot(s_coo_new,psi_curve_new)
% % % % hold on
% % 
% % % % SP_candidates_s = NaN(n_piece,2);
% % % % SP_candidates_psi = NaN(n_piece,2);
% % 
% % SP_candidates_s = [];
% % SP_candidates_psi = [];
% % 
% % vec_n_piece = [10 11];
% % for jj = 1:2
% %     n_piece = vec_n_piece(jj);
% %     s_coo_new = linspace(s_coo(1),s_coo(end),3*n_piece+1).';
% %     psi_curve_new = interp1(s_coo,psi_curve_sort,s_coo_new);
% %     prod_scal_new = interp1(s_coo,prod_scal_sort,s_coo_new);
% %     
% %     % % plot(curve(prod_scal>=0,1),curve(prod_scal>=0,2),'*')
% %     % % quiver(curve(:,1),curve(:,2),vers_n(:,1),vers_n(:,2),'k')
% %     
% %     % % figure
% %     % % plot3(curve(:,1),curve(:,2),psi_curve,'o'); hold on;
% %     % % plot3(curve(prod_scal>0,1),curve(prod_scal>0,2),psi_curve(prod_scal>0),'*'); hold on;
% %     
% %     % % plot(Axis(1),Axis(2),'o'); hold on;
% %     % % plot(curve(prod_scal>0,1),curve(prod_scal>0,2),'*'); hold on;
% %     % % quiver(curve(prod_scal<=0,1),curve(prod_scal<=0,2),grad_curve(prod_scal<=0,1),grad_curve(prod_scal<=0,2),'k'); hold on;
% %     % % quiver(curve(:,1),curve(:,2),grad_curve(:,1),grad_curve(:,2)); hold on;
% %     
% %     
% %     % % theta_new_last = theta_new(end);
% %     % % psi_curve_new_last = psi_curve_new(end);
% %     % %
% %     % % theta_new = theta_new(1:end-1);
% %     % % psi_curve_new = psi_curve_new(1:end-1);
% %     
% %     s_coo_rshp = reshape(s_coo_new(1:end-1),3,numel(s_coo_new(1:end-1))/3).';
% %     prod_scal_rshp = reshape(prod_scal_new(1:end-1),3,numel(prod_scal_new(1:end-1))/3).';
% %     Psi_rshp = reshape(psi_curve_new(1:end-1),3,numel(s_coo_new(1:end-1))/3).';
% %     
% %     s_coo_rshp = [s_coo_rshp [s_coo_rshp(2:end,1); s_coo_new(end)]].';
% %     prod_scal_rshp = [prod_scal_rshp [prod_scal_rshp(2:end,1); prod_scal_new(end)]].';
% %     Psi_rshp = [Psi_rshp [Psi_rshp(2:end,1); psi_curve_new(end)]].';
% %     
% %     
% %     s_coo_rshp = reshape(s_coo_new(1:end-1),3,numel(s_coo_new(1:end-1))/3).';
% %     prod_scal_rshp = reshape(prod_scal_new(1:end-1),3,numel(prod_scal_new(1:end-1))/3).';
% %     Psi_rshp = reshape(psi_curve_new(1:end-1),3,numel(s_coo_new(1:end-1))/3).';
% %     
% %     s_coo_rshp = [s_coo_rshp [s_coo_rshp(2:end,1); s_coo_new(end)]].';
% %     prod_scal_rshp = [prod_scal_rshp [prod_scal_rshp(2:end,1); prod_scal_new(end)]].';
% %     Psi_rshp = [Psi_rshp [Psi_rshp(2:end,1); psi_curve_new(end)]].';
% %     
% %     
% %     
% %     for ii = 1:n_piece
% %         coeffs = fun_spline_1D_coeffs_f_vals(Psi_rshp(:,ii),s_coo_rshp(:,ii));
% %         s_ii = linspace(s_coo_rshp(1,ii),s_coo_rshp(end,ii),20).';
% %         
% %         ft = fun_eval_spline_1D(coeffs,s_ii);
% %         % %     pause
% %         plot(s_ii,ft,'o')
% %         
% %         a = 3*coeffs(4);
% %         b = 2*coeffs(3);
% %         c = coeffs(2);
% %         delta = b^2 - 4*a*c;
% %         x_12 = [-b+sqrt(delta) -b-sqrt(delta)]/2/a;
% %         
% %         x_12(find(imag(x_12))) = NaN;
% %         
% %         x_12(x_12 < s_ii(1)) = NaN;
% %         x_12(x_12 > s_ii(end)) = NaN;
% %         
% %         prod_scal_ii = interp1(s_coo_rshp(:,ii),prod_scal_rshp(:,ii),x_12);
% %         x_12(prod_scal_ii>0) = NaN;
% %         
% %         x_12 = x_12.';
% %         
% %         SP_candidates_s = [SP_candidates_s; x_12];
% %         SP_candidates_psi = [SP_candidates_psi; fun_eval_spline_1D(coeffs,x_12)];
% % % %         plot(x_12,fun_eval_spline_1D(coeffs,x_12),'*r', 'LineWidth',2)
% %         
% %     end
% %     
% % end




%% strategy #2

n_piece = 40;
s_coo_new = linspace(s_coo(1),2*s_coo(end),2*n_piece+1).';

qq1 = [s_coo; s_coo(2:end)+s_coo(end)];
qq2 = [psi_curve_sort; psi_curve_sort(2:end)];
qq3 = [prod_scal_sort; prod_scal_sort(2:end)];
qq4= [curve_sort; curve_sort(2:end,:)];



psi_curve_new = interp1(qq1,qq2,s_coo_new);
prod_scal_new = interp1(qq1,qq3,s_coo_new);

ind_rshp = [ [1 2 3 4];...
    vec_linspace((2:2:2*n_piece-1-1).',(2:2:2*n_piece-1-1).'+3,4)];

s_coo_rshp = zeros(4,n_piece);
Psi_rshp = zeros(4,n_piece);
prod_scal_rshp = zeros(4,n_piece);
for ii = 1:size(ind_rshp,1)
    s_coo_rshp(:,ii) = s_coo_new(ind_rshp(ii,:));
    Psi_rshp(:,ii) = psi_curve_new(ind_rshp(ii,:));
    prod_scal_rshp(:,ii) = prod_scal_new(ind_rshp(ii,:));
end

% % SP_candidates_s = NaN(n_piece,2);
% % SP_candidates_psi = NaN(n_piece,2);


% % figure
% % plot(s_coo_new,psi_curve_new)
% % hold on

SP_candidates_s = [];
SP_candidates_psi = [];

for ii = 1:n_piece
    coeffs = fun_spline_1D_coeffs_f_vals(Psi_rshp(:,ii),s_coo_rshp(:,ii));
    s_ii = linspace(s_coo_rshp(1,ii),s_coo_rshp(end,ii),20).';

    ft = fun_eval_spline_1D(coeffs,s_ii);
    % %     pause
    % %     plot(s_ii,ft,'o')
    
    a = 3*coeffs(4);
    b = 2*coeffs(3);
    c = coeffs(2);
    delta = b^2 - 4*a*c;
    x_12 = [-b+sqrt(delta) -b-sqrt(delta)]/2/a;
    
    x_12(find(imag(x_12))) = NaN;
    
    x_12(x_12 < s_ii(1)) = NaN;
    x_12(x_12 > s_ii(end)) = NaN;
    
    prod_scal_ii = interp1(s_coo_rshp(:,ii),prod_scal_rshp(:,ii),x_12);
    x_12(prod_scal_ii>0) = NaN;
    
    x_12 = x_12.';
    
    SP_candidates_s = [SP_candidates_s; x_12];
    SP_candidates_psi = [SP_candidates_psi; fun_eval_spline_1D(coeffs,x_12)];
% %     plot(x_12,fun_eval_spline_1D(coeffs,x_12),'*r', 'LineWidth',2)
% %     pause
        
end    

s_candidates = SP_candidates_s(:);
s_candidates = s_candidates(~isnan(s_candidates));

psi_candidates = SP_candidates_psi(:);
psi_candidates = psi_candidates(~isnan(psi_candidates));

% %     figure
% %     plot3(curve(:,1),curve(:,2),psi_curve,'ro'); hold on;
% %     plot3(curve(prod_scal<=0,1),curve(prod_scal<=0,2),psi_curve(prod_scal<=0,1),'ok'); hold on;
% %
% %     for ii = 1:numel(psi_candidates)
% %         plot3(interp1(s_coo,curve_sort(:,1),s_candidates(ii)), interp1(s_coo,curve_sort(:,2),s_candidates(ii)), ...
% %             psi_candidates(ii),'*', 'LineWidth',5); hold on;
% %         pause
% %     end

[Psi_Lim,ind_Lim] = max(psi_candidates);
s_lim = s_candidates(ind_Lim);

RZ_Lim = [interp1(qq1,qq4(:,1),s_lim) interp1(qq1,qq4(:,2),s_lim)];

% % rho_lim = interp1(theta_sort,rho_curve_sort,s_candidates(ind_Lim));
% % RZ_Lim = rho_lim*[cos(s_candidates(ind_Lim)) sin(s_candidates(ind_Lim))] + sum(curve)/size(curve,1);



% % plot3(RZ_Lim(1),RZ_Lim(2),Psi_Lim,'*', 'LineWidth',2);
% % 
% % 
% % plot(RZ_Lim(1),RZ_Lim(2),'*', 'LineWidth',2);





end


function OUT = vec_linspace(vec_start,vec_end,nstep)

OUT = zeros(size(vec_end,1),nstep);
for ii=1:size(vec_end,1)
    OUT(ii,:) = linspace(vec_start(ii), vec_end(ii), nstep);
end

end














