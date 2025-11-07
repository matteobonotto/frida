function [solk1] = fun_FluxGauss(meshData,solk1,SETTINGS)

Psi_nodes        = solk1.Psi;
tri              = meshData.t;
shape_functions  = meshData.shape_functions;
N_order          = meshData.shape_functions_N_order;
n_Gauss          = meshData.n_Gauss;
P_Gauss          = meshData.nodes_Gauss_pla;

if SETTINGS.RUN_MEX_ROUTINE == false
    [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
end
solk1.Psi_Gauss = Psi_Gauss;
