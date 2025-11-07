function coeffs = fun_spline_1D_coeffs_f_vals(f_vals,PP)

bb = f_vals;

if size(PP,1) == 1
    PP = PP.';
end

AA = [ones(size(PP,1),1) PP PP.^2 PP.^3];
coeffs = AA\bb;

end



% % bb = [f_vals; fprime_vals*(t_lim(2)-t_lim(1))];
% %
% %
% % coeffs = AA_inv*bb;

