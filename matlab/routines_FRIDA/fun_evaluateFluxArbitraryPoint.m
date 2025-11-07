function [f_rz] = fun_evaluateFluxArbitraryPoint(meshData,f_rz_nodes,Points)

%%

N_order = meshData.shape_functions_N_order;

nodes = meshData.n;
triangles = meshData.t;

rr = meshData.n(:,1);
zz = meshData.n(:,2);
midpoints = [sum(rr(meshData.t),2)/size(meshData.t,2)  sum(zz(meshData.t),2)/size(meshData.t,2)];


npoints = size(Points,1);
ind_tri_sel = zeros(npoints,1);

for ii=1:npoints
    
   P_ii = Points(ii,:);
   
   dist = sqrt((P_ii(1)-midpoints(:,1)).^2 + (P_ii(2)-midpoints(:,2)).^2);
    
   [~,ind_sort] = sort(dist);
   
   INSIDE = 0;
   jj = 0;
   
   while INSIDE == 0
   
       jj = jj+1;
       
       tri_jj = nodes(triangles(ind_sort(jj),1:3),:);
       
       INSIDE = inpolygon(P_ii(1),P_ii(2),tri_jj(:,1),tri_jj(:,2));
       
   end
   
   ind_tri_sel(ii) = ind_sort(jj);
   
    
end
   
    
% % figure
% % triplot(meshData.t(ind_tri_sel,1:3),meshData.n(:,1),meshData.n(:,2),'k'); hold on;
% % plot(Points(:,1),Points(:,2),'o')


%%
tri = meshData.t(ind_tri_sel,:);
shape_functions = meshData.shape_functions(ind_tri_sel,:);


% è possibile sostituirli con i nodi di Gauss
f_nodes = f_rz_nodes(tri);


f_nodes_rep = f_nodes;
% % f_nodes_rep = reshape(f_nodes_rep',3*N_order,3*N_order*size(f_nodes,1))';

f_rz = zeros(size(Points,1),1);


for ii=1:size(tri,2)
    
    if N_order == 1
        coeffs = shape_functions(:,(ii-1)*N_order*3+1:ii*N_order*3);
        
        coeffs_rep = repmat(coeffs,1,3);
% %         coeffs_rep = reshape(coeffs_rep',3,3*size(coeffs,1))';
        
        f_Gauss_ii = f_nodes_rep(:,ii).*(coeffs_rep(:,1).*Points(:,1) + ...
            coeffs_rep(:,2).*Points(:,2) + coeffs_rep(:,3));
        
        f_rz = f_rz+f_Gauss_ii;
        
    elseif N_order == 2
        coeffs = shape_functions(:,(ii-1)*N_order*3+1:ii*N_order*3);
        
        coeffs_rep = repmat(coeffs,1,1);
% %         coeffs_rep = reshape(coeffs_rep',3*N_order,3*N_order*size(coeffs,1))';
        
        fac_geo = coeffs_rep(:,1).*Points(:,1).^2 + ...
            coeffs_rep(:,2).*Points(:,2).^2 + ...
            coeffs_rep(:,3).*Points(:,1).*Points(:,2) + ...
            coeffs_rep(:,4).*Points(:,1) + ...
            coeffs_rep(:,5).*Points(:,2) + ...
            coeffs_rep(:,6);
        
        f_Gauss_ii = f_nodes_rep(:,ii).*fac_geo;
        
        f_rz = f_rz+f_Gauss_ii;

    end
    
end



end