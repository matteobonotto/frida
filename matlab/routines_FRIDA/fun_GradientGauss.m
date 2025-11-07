function [solk1] = fun_GradientGauss(meshData,solk1,SETTINGS)

Psi_nodes        = solk1.Psi;
tri              = meshData.t;
shape_functions  = meshData.shape_functions;
N_order          = meshData.shape_functions_N_order;
n_Gauss          = meshData.n_Gauss;
P_Gauss          = meshData.nodes_Gauss_pla;

if SETTINGS.RUN_MEX_ROUTINE == false
    [Grad_Gauss] = fun_calcGradientGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Grad_Gauss] = fun_calcGradientGauss_FEM_mex(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
end
solk1.Grad_Gauss = Grad_Gauss;