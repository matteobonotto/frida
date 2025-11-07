function [solk1] = fun_GradientNodes(meshData,solk1,SETTINGS)

Psi_nodes        = solk1.Psi;
tri              = meshData.t;
MatInd_nodes_tri = meshData.MatInd_nodes_tri;
nodes            = meshData.n;
shape_functions  = meshData.shape_functions;
N_order          = meshData.shape_functions_N_order;

if SETTINGS.RUN_MEX_ROUTINE == false
    [Grad_nodes] = fun_calcGradientNodes_FEM(tri,nodes,Psi_nodes,MatInd_nodes_tri,shape_functions,N_order);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Grad_nodes] = fun_calcGradientNodes_FEM_mex(tri,nodes,Psi_nodes,MatInd_nodes_tri,shape_functions,N_order);
end

solk1.Grad_nodes = Grad_nodes;