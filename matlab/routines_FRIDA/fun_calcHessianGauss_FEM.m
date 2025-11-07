    function  [Hessian_Psi] = fun_calcHessianGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss)
 
%%
 
nt = size(tri,1);

% % nodes_Gauss_pla_resh = reshape(P_Gauss',2*n_Gauss,nt)';

% % dd_psi_drr_Gauss = zeros(n_Gauss*nt,1);
% % dd_psi_drz_Gauss = zeros(n_Gauss*nt,1);
% % dd_psi_dzz_Gauss = zeros(n_Gauss*nt,1);
   
if N_order == 1
     
     error('second derivatives cannot be used for 1st order basis functions')
     
elseif N_order == 2
     
     
    coeffs_1 = shape_functions(:,1:6);
    coeffs_2 = shape_functions(:,7:12);
    coeffs_3 = shape_functions(:,13:18);
    coeffs_4 = shape_functions(:,19:24);
    coeffs_5 = shape_functions(:,25:30);
    coeffs_6 = shape_functions(:,31:36);
                 
    psi_1 = Psi_nodes(tri(:,1));
    psi_2 = Psi_nodes(tri(:,2));
    psi_3 = Psi_nodes(tri(:,3));
    psi_4 = Psi_nodes(tri(:,4));
    psi_5 = Psi_nodes(tri(:,5));
    psi_6 = Psi_nodes(tri(:,6));
     
    dd_psi_drr_Gauss = zeros(nt,n_Gauss);
    dd_psi_drz_Gauss = zeros(nt,n_Gauss);
    dd_psi_dzz_Gauss = zeros(nt,n_Gauss);
     
    for ii = 1:n_Gauss
                 
        tmp_ii  = 2.*coeffs_1(:,1).*psi_1 + ...
            2.*coeffs_2(:,1).*psi_2 + ...
            2.*coeffs_3(:,1).*psi_3 + ...
            2.*coeffs_4(:,1).*psi_4 + ...
            2.*coeffs_5(:,1).*psi_5 + ...
            2.*coeffs_6(:,1).*psi_6;
        
        dd_psi_drr_Gauss(:,ii) = tmp_ii;
        
        tmp_ii  = coeffs_1(:,3).*psi_1 + ...
            coeffs_2(:,3).*psi_2 + ...
            coeffs_3(:,3).*psi_3 + ...
            coeffs_4(:,3).*psi_4 + ...
            coeffs_5(:,3).*psi_5 + ...
            coeffs_6(:,3).*psi_6;
        
        dd_psi_drz_Gauss(:,ii) = tmp_ii;
        
        tmp_ii  = 2.*coeffs_1(:,2).*psi_1 + ...
            2.*coeffs_2(:,2).*psi_2 + ...
            2.*coeffs_3(:,2).*psi_3 + ...
            2.*coeffs_4(:,2).*psi_4 + ...
            2.*coeffs_5(:,2).*psi_5 + ...
            2.*coeffs_6(:,2).*psi_6;
        
        dd_psi_dzz_Gauss(:,ii) = tmp_ii;

         
    end
     
    dd_psi_drr_Gauss = reshape(dd_psi_drr_Gauss',nt*n_Gauss,1);
    dd_psi_drz_Gauss = reshape(dd_psi_drz_Gauss',nt*n_Gauss,1);
    dd_psi_dzz_Gauss = reshape(dd_psi_dzz_Gauss',nt*n_Gauss,1);
        
% %     quiver(nodes_Gauss_pla(:,1),nodes_Gauss_pla(:,2),grad_psi_r_Gauss,grad_psi_z_Gauss,1)
     
end
 
Hessian_Psi = [dd_psi_drr_Gauss dd_psi_drz_Gauss dd_psi_dzz_Gauss];
 
% % toc
 
 
 
 
 
 
 
 
% % figure
% % plot3(meshData.n(:,1),meshData.n(:,2),Sources,'.')