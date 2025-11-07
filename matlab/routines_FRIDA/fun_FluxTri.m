function [solk1] = fun_FluxTri(meshData,solk1,SETTINGS)

Psi_nodes        = solk1.Psi;
tri              = meshData.t;
shape_functions  = meshData.shape_functions;
N_order          = meshData.shape_functions_N_order;
n_Gauss_fake     = 1;
P_Gauss_fake     = meshData.c_t;

%%%
if SETTINGS.RUN_MEX_ROUTINE == false
    [Psi_c_t] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss_fake,n_Gauss_fake,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Psi_c_t] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss_fake,n_Gauss_fake,shape_functions);
end
solk1.Psi_c_t = Psi_c_t;