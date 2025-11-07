function [ind_in,ind_out] = fun_find_points_INandOUT_separatrix(Points,psi_P,Grad_psi_P,psi_B,Point_Axis)

%% find points inside plasma (psi>psi_b & grad(psi)*vec_Or<0)
% see R. Aòbanese PhD thesis

%%% psi>psi_b
ind_sel_1 = (psi_P >= psi_B);


%%% grad(psi)*vec_Or<0
vec_Or = Points - Point_Axis;
prod_scal = sum(Grad_psi_P.*vec_Or,2);
ind_sel_2 = (prod_scal <= 0);


%%% itersection
ind_in = (ind_sel_1 & ind_sel_2);
ind_out = (true(size(psi_P)) & ~ind_in);

ind_in = find(ind_in);
ind_out = find(ind_out);

end























