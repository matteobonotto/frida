function [s_star,fs_star] = fun_find_SP_1D(s,fs)


%%
if abs((fs(end)-fs(1))/fs(1)) < 1e-5
    s_extended = [s; s(2:3)+s(end)];
    fs_extended = [fs; fs(2:3)];    
else
    s_extended = [s; s(1:2)+s(end)];
    fs_extended = [fs; fs(1:2)];
end


n_piece = 50;
s_new = linspace(s_extended(1),s_extended(end),2*n_piece+1).';
fs_new = interp1(s_extended,fs_extended,s_new);

ind_rshp = [ [1 2 3 4];...
    vec_linspace((2:2:2*n_piece-1-1).',(2:2:2*n_piece-1-1).'+3,4)];

s_coo_rshp = zeros(4,n_piece);
fs_rshp = zeros(4,n_piece);
for ii = 1:size(ind_rshp,1)
    s_coo_rshp(:,ii) = s_new(ind_rshp(ii,:));
    fs_rshp(:,ii) = fs_new(ind_rshp(ii,:));
end

%%

SP_candidates_s = NaN(2*n_piece,1);
SP_candidates_psi = NaN(2*n_piece,1);
ii_pos = 0;

for ii = 1:n_piece
    coeffs = fun_spline_1D_coeffs_f_vals(fs_rshp(:,ii),s_coo_rshp(:,ii));
    s_ii = linspace(s_coo_rshp(1,ii),s_coo_rshp(end,ii),20).';

% %     ft = fun_eval_spline_1D(coeffs,s_ii);
    % %     pause
    % %     plot(s_ii,ft,'o')
    
    a = 3*coeffs(4);
    b = 2*coeffs(3);
    c = coeffs(2);
    delta = b^2 - 4*a*c;
    x_12 = [-b+sqrt(delta) -b-sqrt(delta)]/2/a;
    
    x_12(imag(x_12) ~= 0) = NaN;
    
    x_12(x_12 < s_ii(1)) = NaN;
    x_12(x_12 > s_ii(end)) = NaN;
       
    x_12 = x_12.';
    
    ii_pos = ii_pos(end)+1:ii_pos(end)+length(x_12);

    SP_candidates_s(ii_pos) = x_12;
    SP_candidates_psi(ii_pos) = fun_eval_spline_1D(coeffs,x_12);
% %     plot(x_12,fun_eval_spline_1D(coeffs,x_12),'*r', 'LineWidth',2)
% %     pause
        
end    


%%
ind = ~isnan(SP_candidates_s);

s_star = SP_candidates_s(ind);
fs_star = SP_candidates_psi(ind);

s_star = s_star(s_star <= s(end));
fs_star = fs_star(s_star <= s(end));

























