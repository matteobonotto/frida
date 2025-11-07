function [solk1] = fun_GradientTri(meshData,solk1,SETTINGS)

Psi_nodes        = solk1.Psi;
tri              = meshData.t;
shape_functions  = meshData.shape_functions;
N_order          = meshData.shape_functions_N_order;
n_Gauss_fake     = 1;
P_Gauss_fake     = meshData.c_t;

%%%
if SETTINGS.RUN_MEX_ROUTINE == false
    [Grad_c_t] = fun_calcGradientGauss_FEM(tri,Psi_nodes,P_Gauss_fake,shape_functions,N_order,n_Gauss_fake);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Grad_c_t] = fun_calcGradientGauss_FEM_mex(tri,Psi_nodes,P_Gauss_fake,shape_functions,N_order,n_Gauss_fake);
end
solk1.Grad_c_t = Grad_c_t;