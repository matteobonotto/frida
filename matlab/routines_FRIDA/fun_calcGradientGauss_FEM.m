function  [Gradient_Psi] = fun_calcGradientGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss)
 
%%
 
nt = size(tri,1);

nodes_Gauss_pla_resh = reshape(P_Gauss',2*n_Gauss,nt)';

grad_psi_r_Gauss = zeros(n_Gauss*nt,1);
grad_psi_z_Gauss = zeros(n_Gauss*nt,1);
   
if N_order == 1
     
    coeffs_1 = shape_functions(:,1:3);
    coeffs_2 = shape_functions(:,4:6);
    coeffs_3 = shape_functions(:,7:9);
     
    psi_1 = Psi_nodes(tri(:,1));
    psi_2 = Psi_nodes(tri(:,2));
    psi_3 = Psi_nodes(tri(:,3));
     
    grad_psi_r_Gauss = coeffs_1(:,1).*psi_1 + coeffs_2(:,1).*psi_2 + coeffs_3(:,1).*psi_3;
    grad_psi_z_Gauss = coeffs_1(:,2).*psi_1 + coeffs_2(:,2).*psi_2 + coeffs_3(:,2).*psi_3;
     
elseif N_order == 2
     
     
    coeffs_1 = shape_functions(:,1:6);
    coeffs_2 = shape_functions(:,7:12);
    coeffs_3 = shape_functions(:,13:18);
    coeffs_4 = shape_functions(:,19:24);
    coeffs_5 = shape_functions(:,25:30);
    coeffs_6 = shape_functions(:,31:36);
     
    PP = nodes_Gauss_pla_resh;
             
    psi_1 = Psi_nodes(tri(:,1));
    psi_2 = Psi_nodes(tri(:,2));
    psi_3 = Psi_nodes(tri(:,3));
    psi_4 = Psi_nodes(tri(:,4));
    psi_5 = Psi_nodes(tri(:,5));
    psi_6 = Psi_nodes(tri(:,6));
     
    grad_psi_r_Gauss = zeros(nt,n_Gauss);
    grad_psi_z_Gauss = zeros(nt,n_Gauss);
     
    for ii = 1:n_Gauss
         
        PP_ii = PP(:,(ii-1)*2+1:ii*2);
        tmp_ii  = (2.*coeffs_1(:,1).*PP_ii(:,1) + coeffs_1(:,3).*PP_ii(:,2) + coeffs_1(:,4)).*psi_1 + ...
            (2.*coeffs_2(:,1).*PP_ii(:,1) + coeffs_2(:,3).*PP_ii(:,2) + coeffs_2(:,4)).*psi_2 + ...
            (2.*coeffs_3(:,1).*PP_ii(:,1) + coeffs_3(:,3).*PP_ii(:,2) + coeffs_3(:,4)).*psi_3 + ...
            (2.*coeffs_4(:,1).*PP_ii(:,1) + coeffs_4(:,3).*PP_ii(:,2) + coeffs_4(:,4)).*psi_4 + ...
            (2.*coeffs_5(:,1).*PP_ii(:,1) + coeffs_5(:,3).*PP_ii(:,2) + coeffs_5(:,4)).*psi_5 + ...
            (2.*coeffs_6(:,1).*PP_ii(:,1) + coeffs_6(:,3).*PP_ii(:,2) + coeffs_6(:,4)).*psi_6;
         
        grad_psi_r_Gauss(:,ii) = tmp_ii;
         
        tmp_ii  = (2.*coeffs_1(:,2).*PP_ii(:,2) + coeffs_1(:,3).*PP_ii(:,1) + coeffs_1(:,5)).*psi_1 + ...
            (2.*coeffs_2(:,2).*PP_ii(:,2) + coeffs_2(:,3).*PP_ii(:,1) + coeffs_2(:,5)).*psi_2 + ...
            (2.*coeffs_3(:,2).*PP_ii(:,2) + coeffs_3(:,3).*PP_ii(:,1) + coeffs_3(:,5)).*psi_3 + ...
            (2.*coeffs_4(:,2).*PP_ii(:,2) + coeffs_4(:,3).*PP_ii(:,1) + coeffs_4(:,5)).*psi_4 + ...
            (2.*coeffs_5(:,2).*PP_ii(:,2) + coeffs_5(:,3).*PP_ii(:,1) + coeffs_5(:,5)).*psi_5 + ...
            (2.*coeffs_6(:,2).*PP_ii(:,2) + coeffs_6(:,3).*PP_ii(:,1) + coeffs_6(:,5)).*psi_6;
         
        grad_psi_z_Gauss(:,ii) = tmp_ii;
         
    end
     
    grad_psi_r_Gauss = reshape(grad_psi_r_Gauss',nt*n_Gauss,1);
    grad_psi_z_Gauss = reshape(grad_psi_z_Gauss',nt*n_Gauss,1);
        
% %     quiver(nodes_Gauss_pla(:,1),nodes_Gauss_pla(:,2),grad_psi_r_Gauss,grad_psi_z_Gauss,1)
     
end
 
Gradient_Psi = [grad_psi_r_Gauss grad_psi_z_Gauss];
 
% % toc
 
 
 
 
 
 
 
 
% % figure
% % plot3(meshData.n(:,1),meshData.n(:,2),Sources,'.')