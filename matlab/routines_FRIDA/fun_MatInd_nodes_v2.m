function [MatInd_nodes] = fun_MatInd_nodes_v2(tri,nodes,MatInd_nodes_tri)


%% Incidence Matrix between triangles and Gauss points

nn = size(nodes,1);

MatInd_nodes = zeros(nn,20);

for ii = 1:nn
    
    matind_ii = MatInd_nodes_tri(ii,MatInd_nodes_tri(ii,:) ~= 0);
    
    tri_vicini = tri(matind_ii,:);
    
    matind_nodes = unique(tri_vicini(:));
    
    MatInd_nodes(ii,1:numel(matind_nodes)) = matind_nodes;
    
end


