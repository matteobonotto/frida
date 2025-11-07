function [KK_nobc] = fun_StiffMat_FEM(meshData_loc,MEX_OPT)

tri              = meshData_loc.t;
nodes            = meshData_loc.n;
N_order          = meshData_loc.shape_functions_N_order;
GAUSS_QUAD_DEGREE_STIFFMAT          = meshData_loc.GAUSS_QUAD_DEGREE_STIFFMAT;
shape_functions  = meshData_loc.shape_functions;

if MEX_OPT == false
    [irow,jcol,KK_vals] = fun_StiffMat_FEM_fast(tri,nodes,N_order,shape_functions,GAUSS_QUAD_DEGREE_STIFFMAT);
    
elseif MEX_OPT == true
    [irow,jcol,KK_vals] = fun_StiffMat_FEM_fast_mex(tri,nodes,N_order,shape_functions,GAUSS_QUAD_DEGREE_STIFFMAT);
    
end

KK_nobc = sparse(irow,jcol,KK_vals);