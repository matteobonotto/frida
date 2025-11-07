function ft = fun_eval_spline_1D(coeffs,t)

ft = coeffs(1,:) + coeffs(2,:).*t + coeffs(3,:).*t.^2 + coeffs(4,:).*t.^3;
