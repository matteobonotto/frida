function [Grad_nodes] = fun_calcGradientNodes_FEM(tri,nodes,Psi_nodes,MatInd_nodes_tri,shape_functions,N_order)

nn = length(Psi_nodes);

Grad_r = zeros(nn,1);
Grad_z = zeros(nn,1);

for ii=1:nn % parfor nella versione mexata
    
   ind_t_ii = MatInd_nodes_tri(ii,1);

   tri_ii = tri(ind_t_ii,:);
   
   r_node = nodes(ii,1);
   z_node = nodes(ii,2);
   Psi_node = Psi_nodes(tri_ii);
   
   shape_functions_ii = shape_functions(ind_t_ii,:);
   
   cc = reshape(shape_functions_ii,3*N_order,3*N_order)';
      
   Grad_r_ii = 0;
   Grad_z_ii = 0;
   
   for jj = 1:3*N_order
       
       Grad_r_jj = (2*cc(jj,1)*r_node + cc(jj,3)*z_node + cc(jj,4))*Psi_node(jj);
       Grad_z_jj = (2*cc(jj,2)*z_node + cc(jj,3)*r_node + cc(jj,5))*Psi_node(jj);
       
       Grad_r_ii = Grad_r_ii + Grad_r_jj;
       Grad_z_ii = Grad_z_ii + Grad_z_jj;
       
   end
   
   Grad_r(ii) = Grad_r_ii;
   Grad_z(ii) = Grad_z_ii;
      
end


Grad_nodes = [Grad_r Grad_z];