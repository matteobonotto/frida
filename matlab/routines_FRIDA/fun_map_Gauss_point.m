function [P_G_mapped,w_mapped] = fun_map_Gauss_point(nodes_norm,ww_norm,n_Gauss,P1,P2,P3)


% mapping from unit triangle [(0,0),(1,0),(0,1)] to arbitrary triangle [P1,P2,P3]
T = [P2' P3'] - [P1' P1'];

tmp1 = repmat(P1',1,n_Gauss);
tmp2 = T*nodes_norm';
nodes_ii =  (tmp1 + tmp2);
P_G_mapped =  nodes_ii';

% weigths (ww_norm is normalized by triangle area)
edge_1 = P2 - P1;
edge_2 = P3 - P2;

% %         polyarea([P1(1) P2(1) P3(1)],[P1(2) P2(2) P3(2)])
det_Jac = .5*abs(edge_1(1)*edge_2(2) - edge_1(2)*edge_2(1));

% %     area = polyarea(triangle(:,1),triangle(:,2));
area = det_Jac;
w_mapped = area*ww_norm;
