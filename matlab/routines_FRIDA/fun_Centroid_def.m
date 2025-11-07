function [r_C,z_C,Psi_C] = fun_Centroid_def(meshData_loc,solk)


chi_0 = 1;
q_0 = (chi_0.*solk.Jphi).'*meshData_loc.ww_Gauss_pla;


chi_1 = meshData_loc.nodes_Gauss_pla(:,2);
q_1 = (chi_1.*solk.Jphi).'*meshData_loc.ww_Gauss_pla;
z_C = q_1/q_0;


chi_2 = (meshData_loc.nodes_Gauss_pla(:,1)).^2;
q_2 = (chi_2.*solk.Jphi).'*meshData_loc.ww_Gauss_pla;
r_C = sqrt(q_2/q_0);

%%

distanza=sqrt((r_C - meshData_loc.nodes_Gauss_pla(:,1)).^2+...
    (z_C - meshData_loc.nodes_Gauss_pla(:,2)).^2);
[~,ind_min]=min(distanza);

ind_t_Centroid = ceil(ind_min/meshData_loc.n_Gauss);
ind_n_centroid = meshData_loc.t(ind_t_Centroid,:);

N_order = meshData_loc.shape_functions_N_order;
cc_tri = reshape(meshData_loc.shape_functions(ind_t_Centroid,:),3*N_order,3*N_order).';

Psi_C = solk.Psi(ind_n_centroid).'*(cc_tri*([r_C^2 z_C^2 r_C*z_C r_C z_C 1].'));


