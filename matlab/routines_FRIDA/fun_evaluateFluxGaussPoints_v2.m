function [f_rz_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,f_rz_nodes,N_order,P_Gauss,n_Gauss,shape_functions)


nt = size(tri,1);

f_nodes = f_rz_nodes(tri);


f_nodes_rep = repmat(f_nodes,1,n_Gauss);
f_nodes_rep = reshape(f_nodes_rep',3*N_order,n_Gauss*nt)';

f_rz_Gauss_temp = zeros(size(P_Gauss,1),N_order*3);


for ii=1:size(tri,2) %parfor nella versione mexata
    
    if N_order == 1
        coeffs = shape_functions(:,(ii-1)*N_order*3+1:ii*N_order*3);
        
        coeffs_rep = repmat(coeffs,1,n_Gauss);
        coeffs_rep = reshape(coeffs_rep',3*N_order,n_Gauss*size(coeffs,1))';
        
        f_Gauss_ii = f_nodes_rep(:,ii).*(coeffs_rep(:,1).*P_Gauss(:,1) + ...
            coeffs_rep(:,2).*P_Gauss(:,2) + coeffs_rep(:,3));
        
        f_rz_Gauss_temp(:,ii) = f_Gauss_ii;
        
    elseif N_order == 2
        coeffs = shape_functions(:,(ii-1)*N_order*3+1:ii*N_order*3);
        
        coeffs_rep = repmat(coeffs,1,n_Gauss);
        coeffs_rep = reshape(coeffs_rep',3*N_order,n_Gauss*size(coeffs,1))';
        
        fac_geo = coeffs_rep(:,1).*P_Gauss(:,1).^2 + ...
            coeffs_rep(:,2).*P_Gauss(:,2).^2 + ...
            coeffs_rep(:,3).*P_Gauss(:,1).*P_Gauss(:,2) + ...
            coeffs_rep(:,4).*P_Gauss(:,1) + ...
            coeffs_rep(:,5).*P_Gauss(:,2) + ...
            coeffs_rep(:,6);
        
        f_Gauss_ii = f_nodes_rep(:,ii).*fac_geo;
        
        f_rz_Gauss_temp(:,ii) = f_Gauss_ii;

    end
    
end

%
f_rz_Gauss = sum(f_rz_Gauss_temp,2);









