function res_int = fun_integrate_poly_on_tri(exp_r,exp_z)

res_int = factorial(exp_r)*factorial(exp_z)/factorial(2+exp_r+exp_z);
% % for jj = 1:length(exp_r)
% %     res_int = res_int + factorial(exp_r(jj))*factorial(exp_z(jj))/factorial(2+exp_r(jj)+exp_z(jj));
% % end