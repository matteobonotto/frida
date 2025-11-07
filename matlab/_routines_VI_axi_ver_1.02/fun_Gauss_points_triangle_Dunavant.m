function [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,GAUSS_QUAD_DEGREE)  %#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GAUSS_QUAD_DEGREE = order of the quadrature rule
% n_Gauss = number of points

% see tab. A. Sommariva https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/rules_triangle.html
% -> DUNAVANT # [Corrected] 
% see D. A. DUNAVANT, "HIGH DEGREE EFFICIENT SYMMETRICAL GAUSSIAN QUADRATURE RULES FOR THE TRIANGLE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Dunavant Gauss points for unit triangle
[xyw] = set_dunavant_C_standard(GAUSS_QUAD_DEGREE);

n_Gauss = size(xyw,1);

nodes_norm = xyw(:,[1 2]);
ww_norm = xyw(:,3);
ww_norm = ww_norm/.5;

%% Map Gauss points into real triangles

nodes = zeros(n_Gauss*size(P1,1),2);
ww = zeros(n_Gauss*size(P1,1),1);

for ii = 1:size(P1,1)
    
    P1_ii = P1(ii,:);
    P2_ii = P2(ii,:);
    P3_ii = P3(ii,:);

    % mapping from unit triangle [(0,0),(1,0),(0,1)] to arbitrary triangle [P1,P2,P3]
    [P_G_mapped,w_mapped] = fun_map_Gauss_point(nodes_norm,ww_norm,n_Gauss,P1_ii,P2_ii,P3_ii);
    
    index_nodes = (ii-1)*n_Gauss+1:ii*n_Gauss;
    
    nodes(index_nodes,1) = P_G_mapped(:,1);
    nodes(index_nodes,2) = P_G_mapped(:,2);
    
    ww(index_nodes) = w_mapped;
    
    
    % %     triplot([1 2 3], [P1_ii(1); P2_ii(1); P3_ii(1)],[P1_ii(2); P2_ii(2); P3_ii(2)])
    % %     plot(PP(1),PP(2),'or'); hold on; axis equal;
end
































