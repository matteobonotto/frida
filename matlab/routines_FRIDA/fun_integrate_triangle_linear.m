function [I_tri,PP_Gauss] = fun_integrate_triangle_linear(...
    F_ref,tri_local,shape_f_jj,nodes_jj,shape_f_norm,n_Gauss_tri,N_order)


%% Integrating on actual element

% % Map Gauss points onto the real triangle 
P1 = tri_local(1,:);
P2 = tri_local(2,:);
P3 = tri_local(3,:);

n_Gauss = n_Gauss_tri;
[ww_G,PP_Gauss] = fun_GaussQuadTri(P1,P2,P3,n_Gauss);

% Integrate
gg_Gauss_jj = F_ref(PP_Gauss(:,1),PP_Gauss(:,2));

tmp_ii = ww_G.*(gg_Gauss_jj);

tmp = [PP_Gauss.^2 PP_Gauss(:,1).*PP_Gauss(:,2) ...
    PP_Gauss ones(n_Gauss_tri,1)]';

cc = reshape(shape_f_jj,3*N_order,3*N_order)';

NN_1_norm = cc(1,:)*tmp;
NN_2_norm = cc(2,:)*tmp;
NN_3_norm = cc(3,:)*tmp;
NN_4_norm = cc(4,:)*tmp;
NN_5_norm = cc(5,:)*tmp;
NN_6_norm = cc(6,:)*tmp;

I_tri_1 = tmp_ii'*NN_1_norm';
I_tri_2 = tmp_ii'*NN_2_norm';
I_tri_3 = tmp_ii'*NN_3_norm';
I_tri_4 = tmp_ii'*NN_4_norm';
I_tri_5 = tmp_ii'*NN_5_norm';
I_tri_6 = tmp_ii'*NN_6_norm';

I_tri = [I_tri_1; I_tri_2; I_tri_3; I_tri_4; I_tri_5; I_tri_6];


%% Integrating on master element

% % % % figure
% % % % plot([0 1 0 0], [0 0 1 0], 'k'); hold on; axis equal
% % % % plot(tri_local_norm([1:end 1],1),tri_local_norm([1:end 1],2),'ro')
% % 
% % % Map Gauss points onto the real quadrilateral (16 points) element
% % PP_ii = nodes_jj(1:3,:);
% % 
% % P1_ii = PP_ii(1,:);
% % P2_ii = PP_ii(2,:);
% % P3_ii = PP_ii(3,:);
% % 
% % T = [P2_ii' P3_ii'] - [P1_ii' P1_ii'];
% % 
% % tri_local_norm = (inv(T)*(tri_local'-P1_ii'))';
% % 
% % % % figure
% % % % plot([0 1 0 0], [0 0 1 0], 'k'); hold on; axis equal
% % % % plot(tri_local_norm([1:end 1],1),tri_local_norm([1:end 1],2),'r')
% % 
% % Jac = abs(T(1,1)*T(2,2) - T(1,2)*T(2,1));
% % 
% % 
% % P1 = tri_local_norm(1,:);
% % P2 = tri_local_norm(2,:);
% % P3 = tri_local_norm(3,:);
% % 
% % n_Gauss = n_Gauss_tri;
% % [ww_G_norm,PP_Gauss_norm] = fun_GaussQuadTri(P1,P2,P3,n_Gauss);
% % 
% % PP_Gauss = (P1_ii' + T*PP_Gauss_norm')';
% % 
% % % Integrate
% % gg_Gauss_jj = F_ref(PP_Gauss(:,1),PP_Gauss(:,2));
% % 
% % tmp_ii = ww_G_norm.*(gg_Gauss_jj)*Jac;
% % 
% % tmp = [PP_Gauss_norm.^2 PP_Gauss_norm(:,1).*PP_Gauss_norm(:,2) ...
% %     PP_Gauss_norm ones(n_Gauss_tri,1)]';
% % 
% % cc = reshape(shape_f_norm,3*N_order,3*N_order)';
% % 
% % NN_1_norm = cc(1,:)*tmp;
% % NN_2_norm = cc(2,:)*tmp;
% % NN_3_norm = cc(3,:)*tmp;
% % NN_4_norm = cc(4,:)*tmp;
% % NN_5_norm = cc(5,:)*tmp;
% % NN_6_norm = cc(6,:)*tmp;
% % 
% % I_tri_1 = tmp_ii'*NN_1_norm';
% % I_tri_2 = tmp_ii'*NN_2_norm';
% % I_tri_3 = tmp_ii'*NN_3_norm';
% % I_tri_4 = tmp_ii'*NN_4_norm';
% % I_tri_5 = tmp_ii'*NN_5_norm';
% % I_tri_6 = tmp_ii'*NN_6_norm';
% % 
% % I_tri = [I_tri_1; I_tri_2; I_tri_3; I_tri_4; I_tri_5; I_tri_6];











%%

% % 
% % figure
% % hold on; axis equal
% % plot(nodes_jj([1:3 1],1),nodes_jj([1:3 1],2),'r-o')
% % plot(tri_local([1:end 1],1),tri_local([1:end 1],2),'k')
% % 
% % P1_ii = tri_local(1,:);
% % P2_ii = tri_local(2,:);
% % P3_ii = tri_local(3,:);
% % 
% % 
% % m_12 = (P2_ii(2)-P1_ii(2))/(P2_ii(1)-P1_ii(1))
% % 
% % if m_12 < eps
% %     P_star = [P1_ii(1) P3_ii(2)];
% %     
% % elseif isinf(m_12)
% %     P_star = [P3_ii(1) P1_ii(2)];
% %     
% % else
% %     
% %     x_star = 1/(m_12 + 1/m_12)*(P3_ii(2)-P2_ii(2)+m_12*P2_ii(1)+1/m_12*P3_ii(1))
% %     y_star = P2_ii(2)+m_12*(x_star-P2_ii(1))
% %     
% %     P_star = [x_star y_star];
% % end
% % 
% % plot(P_star(1),P_star(2),'o')
% % 
% % aa = norm(P_star - P2_ii)
% % bb = norm(P_star - P1_ii)
% % cc = norm(P_star - P3_ii)
% % 
% % 
% % 
% % P1 = [-bb 0];
% % P2 = [aa 0];
% % P3 = [0 cc];
% % P4 = .5*(P1+P2);
% % P5 = .5*(P2+P3);
% % P6 = .5*(P3+P1);
% % 
% % temp_tri = [P1; P2; P3; P4; P5; P6];
% % 
% % figure
% % plot(temp_tri([1:6 1],1),temp_tri([1:6 1],2),'o')
% % 
% % AA = [P1.^2 P1(1)*P1(2) P1 1; ...
% %     P2.^2 P2(1)*P2(2) P2 1; ...
% %     P3.^2 P3(1)*P3(2) P3 1; ...
% %     P4.^2 P4(1)*P4(2) P4 1; ...
% %     P5.^2 P5(1)*P5(2) P5 1; ...
% %     P6.^2 P6(1)*P6(2) P6 1];
% % coeff = inv(AA);
% % 
% % exponents = [2 0; 0 2; 1 1; 1 0; 0 1; 0 0];
% % 
% % m = exponents(:,1)
% % n = exponents(:,2)
% % 
% % tmp_Int = cc.^(n+1).*(aa.^(m+1)-(-bb).^(m+1)).*factorial(m).*factorial(m)./factorial(2+m+n)
% % 
% % coeffs = reshape(shape_f_jj,3*N_order,3*N_order)';
% % 
% % 
% % gg_nodes_jj = F_ref(nodes_jj(:,1),nodes_jj(:,2));
% % coeffs = reshape(shape_f_jj,3*N_order,3*N_order)';
% % gg_nodes_jj(gg_nodes_jj<0) = 0;
% % 
% % % % (coeff*tmp_Int).*gg_nodes_jj
% % 
% % 
% % 2*coeff(4,:)*gg_nodes_jj(4)*tmp_Int
% % 
% % polyarea(tri_local(:,1),tri_local(:,2))
% % 
% % % map subtriangle into master element
% % PP_ii = nodes_jj(1:3,:);
% % 
% % P1_ii = tri_local(1,:);
% % P2_ii = tri_local(2,:);
% % P3_ii = tri_local(3,:);
% % 
% % T = [P2_ii' P3_ii'] - [P1_ii' P1_ii'];
% % 
% % 
% % % Original triangle (mesh element) mapped to have the subtriangle as master
% % % element
% % nodes_jj_local_norm = (inv(T)*(nodes_jj'-P1_ii'))';
% % 
% % 
% % % Shape function of the mapped triangle 
% % 
% % P1 = nodes_jj(1,:);
% % P2 = nodes_jj(2,:);
% % P3 = nodes_jj(3,:);
% % P4 = nodes_jj(4,:);
% % P5 = nodes_jj(5,:);
% % P6 = nodes_jj(6,:);
% % 
% % AA = [P1.^2 P1(1)*P1(2) P1 1; ...
% %     P2.^2 P2(1)*P2(2) P2 1; ...
% %     P3.^2 P3(1)*P3(2) P3 1; ...
% %     P4.^2 P4(1)*P4(2) P4 1; ...
% %     P5.^2 P5(1)*P5(2) P5 1; ...
% %     P6.^2 P6(1)*P6(2) P6 1];
% % coeff = inv(AA);
% % 
% % 
% % % Integrate this shape functions on the master subelement
% % 
% % 
% % 
% % 
% % 
% % figure
% % plot([0 1 0 0], [0 0 1 0], 'k'); hold on; axis equal
% % plot(nodes_jj_local_norm([1:end 1],1),nodes_jj_local_norm([1:end 1],2),'ro')
% % 
% % 
% % 
% % 
% % % Map Gauss points onto the real quadrilateral (16 points) element
% % PP_ii = nodes_jj(1:3,:);
% % 
% % P1_ii = PP_ii(1,:);
% % P2_ii = PP_ii(2,:);
% % P3_ii = PP_ii(3,:);
% % 
% % 
% % T = [P2_ii' P3_ii'] - [P1_ii' P1_ii'];
% % 
% % tri_local_norm = (inv(T)*(tri_local'-P1_ii'))';
% % 
% % figure
% % plot([0 1 0 0], [0 0 1 0], 'k'); hold on; axis equal
% % plot(tri_local_norm([1:end 1],1),tri_local_norm([1:end 1],2),'r')
% % 
% % Jac = abs(T(1,1)*T(2,2) - T(1,2)*T(2,1));
% % 
% % 
% % P1 = tri_local_norm(1,:);
% % P2 = tri_local_norm(2,:);
% % P3 = tri_local_norm(3,:);
% % 
% % n_Gauss = n_Gauss_tri;
% % [ww_G_norm,PP_Gauss_norm] = fun_GaussQuadTri(P1,P2,P3,n_Gauss);
% % 
% % PP_Gauss = (P1_ii' + T*PP_Gauss_norm')';
% % 
% % % Integrate
% % gg_Gauss_jj = F_ref(PP_Gauss(:,1),PP_Gauss(:,2));
% % 
% % tmp_ii = ww_G_norm.*(gg_Gauss_jj)*Jac;
% % 
% % tmp = [PP_Gauss_norm.^2 PP_Gauss_norm(:,1).*PP_Gauss_norm(:,2) ...
% %     PP_Gauss_norm ones(n_Gauss_tri,1)]';
% % 
% % cc = reshape(shape_f_norm,3*N_order,3*N_order)';
% % 
% % NN_1_norm = cc(1,:)*tmp;
% % NN_2_norm = cc(2,:)*tmp;
% % NN_3_norm = cc(3,:)*tmp;
% % NN_4_norm = cc(4,:)*tmp;
% % NN_5_norm = cc(5,:)*tmp;
% % NN_6_norm = cc(6,:)*tmp;
% % 
% % I_tri_1 = tmp_ii'*NN_1_norm';
% % I_tri_2 = tmp_ii'*NN_2_norm';
% % I_tri_3 = tmp_ii'*NN_3_norm';
% % I_tri_4 = tmp_ii'*NN_4_norm';
% % I_tri_5 = tmp_ii'*NN_5_norm';
% % I_tri_6 = tmp_ii'*NN_6_norm';
% % 
% % I_tri = [I_tri_1; I_tri_2; I_tri_3; I_tri_4; I_tri_5; I_tri_6]';




