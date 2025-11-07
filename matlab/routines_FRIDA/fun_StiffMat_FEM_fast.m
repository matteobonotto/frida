function [irow,jcol,KK_vals] = fun_StiffMat_FEM_fast(tri,nodes,N_order,shape_functions,GAUSS_QUAD_DEGREE_STIFFMAT) %#codegen

nt = size(tri,1);

P1 = [0 0];
P2 = [1 0];
P3 = [0 1];

[ww_unit,nodes_G_unit,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,GAUSS_QUAD_DEGREE_STIFFMAT);
ww_unit = ww_unit/.5;

%%
irow    = [];
jcol    = [];
KK_vals = [];

if N_order == 1
        
    irow = zeros((3*N_order)^2*nt,1);        % preallocation: row indices
    jcol = zeros((3*N_order)^2*nt,1);        % preallocation: column indices
    KK_vals = zeros((3*N_order)^2*nt,1);
    
    for ii = 1:nt
        
        tri_ii = tri(ii,:);

        % Map Gauss points from unit triangle [(0,0),(1,0),(0,1)] to arbitrary triangle [P1,P2,P3]
        P_ii = nodes(tri_ii(1:3),:);
        
        P1_ii = P_ii(1,:);
        P2_ii = P_ii(2,:);
        P3_ii = P_ii(3,:);
        
        % mapping 
        T = [P2_ii' P3_ii'] - [P1_ii' P1_ii'];
        
        tmp1 = repmat(P1_ii',1,n_Gauss);
        tmp2 = T*nodes_G_unit';
        tmp3 =  (tmp1 + tmp2);
        P_Gauss_ii =  tmp3';
        
        % weigths (ww_norm is normalized by triangle area)        
        edge_1 = P2_ii - P1_ii;
        edge_2 = P3_ii - P2_ii;
        
        det_Jac = .5*abs(edge_1(1)*edge_2(2) - edge_1(2)*edge_2(1));
        
        area = det_Jac;
        ww_Gauss_ii = area*ww_unit;
        
        
        % Compute integral
        cc = shape_functions(ii,:);
        cc = reshape(cc,3,3*N_order)';
        
        rr_G_ii = P_Gauss_ii(:,1);
        
        fac_Gauss = ww_Gauss_ii'*(1./rr_G_ii);
        
        fac_11 = cc(1,1)*cc(1,1) + cc(1,2)*cc(1,2);
        fac_22 = cc(2,1)*cc(2,1) + cc(2,2)*cc(2,2);
        fac_33 = cc(3,1)*cc(3,1) + cc(3,2)*cc(3,2);
        
        fac_12 = cc(1,1)*cc(2,1) + cc(1,2)*cc(2,2);
        fac_23 = cc(2,1)*cc(3,1) + cc(2,2)*cc(3,2);
        fac_31 = cc(3,1)*cc(1,1) + cc(3,2)*cc(1,2);
                
        kk_11 = fac_11*fac_Gauss;
        kk_22 = fac_22*fac_Gauss;
        kk_33 = fac_33*fac_Gauss;
        
        kk_12 = fac_12*fac_Gauss;
        kk_23 = fac_23*fac_Gauss;
        kk_31 = fac_31*fac_Gauss;
        
        kk_loc = [kk_11 kk_12 kk_31; ...
            kk_12 kk_22 kk_23; ...
            kk_31 kk_23 kk_33];
                
        % Sparse indexing
        num_el = numel(kk_loc);
        idx = num_el*(ii-1)+1:num_el*ii;    % block of 9 indices belonging to i-th triangle
        irow(idx) = repmat(tri_ii',3,1);    % row pointers
        jcol(idx) = reshape(repmat(tri_ii,3,1),9,1);    % column pointers
        
        KK_vals(idx) = kk_loc(:);    % nonzeros
        
    end
    
    
elseif N_order == 2
    
    
    irow = zeros((3*N_order)^2*nt,1);        % preallocation: row indices
    jcol = zeros((3*N_order)^2*nt,1);        % preallocation: column indices
    KK_vals = zeros((3*N_order)^2*nt,1);
    
    for ii = 1:nt
        
        tri_ii = tri(ii,:);
        
        % Map Gauss points from unit triangle [(0,0),(1,0),(0,1)] to arbitrary triangle [P1,P2,P3]
        P_ii = nodes(tri_ii(1:3),:);
        
        P1_ii = P_ii(1,:);
        P2_ii = P_ii(2,:);
        P3_ii = P_ii(3,:);
        
        % mapping 
        T = [P2_ii' P3_ii'] - [P1_ii' P1_ii'];
        
        tmp1 = repmat(P1_ii',1,n_Gauss);
        tmp2 = T*nodes_G_unit';
        tmp3 =  (tmp1 + tmp2);
        P_Gauss_ii =  tmp3';
        
        % weigths (ww_norm is normalized by triangle area)        
        edge_1 = P2_ii - P1_ii;
        edge_2 = P3_ii - P2_ii;
        
        det_Jac = .5*abs(edge_1(1)*edge_2(2) - edge_1(2)*edge_2(1));
        
        area = det_Jac;
        ww_Gauss_ii = area*ww_unit;

        
        % Compute integral
        coeffs = shape_functions(ii,:);
        coeffs = reshape(coeffs,3*N_order,3*N_order)';
        
        aa = coeffs(:,1);
        bb = coeffs(:,2);
        cc = coeffs(:,3);
        dd = coeffs(:,4);
        ee = coeffs(:,5);               
        
        rr_G_ii = P_Gauss_ii(:,1);
        zz_G_ii = P_Gauss_ii(:,2);
                
        gradN_1 = [2*aa(1)*rr_G_ii+cc(1)*zz_G_ii+dd(1)  2*bb(1)*zz_G_ii+cc(1)*rr_G_ii+ee(1)];
        gradN_2 = [2*aa(2)*rr_G_ii+cc(2)*zz_G_ii+dd(2)  2*bb(2)*zz_G_ii+cc(2)*rr_G_ii+ee(2)];
        gradN_3 = [2*aa(3)*rr_G_ii+cc(3)*zz_G_ii+dd(3)  2*bb(3)*zz_G_ii+cc(3)*rr_G_ii+ee(3)];
        gradN_4 = [2*aa(4)*rr_G_ii+cc(4)*zz_G_ii+dd(4)  2*bb(4)*zz_G_ii+cc(4)*rr_G_ii+ee(4)];
        gradN_5 = [2*aa(5)*rr_G_ii+cc(5)*zz_G_ii+dd(5)  2*bb(5)*zz_G_ii+cc(5)*rr_G_ii+ee(5)];
        gradN_6 = [2*aa(6)*rr_G_ii+cc(6)*zz_G_ii+dd(6)  2*bb(6)*zz_G_ii+cc(6)*rr_G_ii+ee(6)];
        
        fac_Gauss = ww_Gauss_ii.*(1./rr_G_ii);
        
        kk_11 = fac_Gauss'*(gradN_1(:,1).*gradN_1(:,1)+gradN_1(:,2).*gradN_1(:,2));
        kk_22 = fac_Gauss'*(gradN_2(:,1).*gradN_2(:,1)+gradN_2(:,2).*gradN_2(:,2));
        kk_33 = fac_Gauss'*(gradN_3(:,1).*gradN_3(:,1)+gradN_3(:,2).*gradN_3(:,2));
        kk_44 = fac_Gauss'*(gradN_4(:,1).*gradN_4(:,1)+gradN_4(:,2).*gradN_4(:,2));
        kk_55 = fac_Gauss'*(gradN_5(:,1).*gradN_5(:,1)+gradN_5(:,2).*gradN_5(:,2));
        kk_66 = fac_Gauss'*(gradN_6(:,1).*gradN_6(:,1)+gradN_6(:,2).*gradN_6(:,2));

        kk_12 = fac_Gauss'*(gradN_1(:,1).*gradN_2(:,1)+gradN_1(:,2).*gradN_2(:,2));
        kk_13 = fac_Gauss'*(gradN_1(:,1).*gradN_3(:,1)+gradN_1(:,2).*gradN_3(:,2));
        kk_14 = fac_Gauss'*(gradN_1(:,1).*gradN_4(:,1)+gradN_1(:,2).*gradN_4(:,2));
        kk_15 = fac_Gauss'*(gradN_1(:,1).*gradN_5(:,1)+gradN_1(:,2).*gradN_5(:,2));
        kk_16 = fac_Gauss'*(gradN_1(:,1).*gradN_6(:,1)+gradN_1(:,2).*gradN_6(:,2));

        kk_23 = fac_Gauss'*(gradN_2(:,1).*gradN_3(:,1)+gradN_2(:,2).*gradN_3(:,2));
        kk_24 = fac_Gauss'*(gradN_2(:,1).*gradN_4(:,1)+gradN_2(:,2).*gradN_4(:,2));
        kk_25 = fac_Gauss'*(gradN_2(:,1).*gradN_5(:,1)+gradN_2(:,2).*gradN_5(:,2));
        kk_26 = fac_Gauss'*(gradN_2(:,1).*gradN_6(:,1)+gradN_2(:,2).*gradN_6(:,2));

        kk_34 = fac_Gauss'*(gradN_3(:,1).*gradN_4(:,1)+gradN_3(:,2).*gradN_4(:,2));
        kk_35 = fac_Gauss'*(gradN_3(:,1).*gradN_5(:,1)+gradN_3(:,2).*gradN_5(:,2));
        kk_36 = fac_Gauss'*(gradN_3(:,1).*gradN_6(:,1)+gradN_3(:,2).*gradN_6(:,2));

        kk_45 = fac_Gauss'*(gradN_4(:,1).*gradN_5(:,1)+gradN_4(:,2).*gradN_5(:,2));
        kk_46 = fac_Gauss'*(gradN_4(:,1).*gradN_6(:,1)+gradN_4(:,2).*gradN_6(:,2));

        kk_56 = fac_Gauss'*(gradN_5(:,1).*gradN_6(:,1)+gradN_5(:,2).*gradN_6(:,2));
        
% %         kk_22 = fac_Gauss'*diag(gradN_2*gradN_2');
% %         kk_33 = fac_Gauss'*diag(gradN_3*gradN_3');
% %         kk_44 = fac_Gauss'*diag(gradN_4*gradN_4');
% %         kk_55 = fac_Gauss'*diag(gradN_5*gradN_5');
% %         kk_66 = fac_Gauss'*diag(gradN_6*gradN_6');
        
% %         kk_12 = fac_Gauss'*diag(gradN_1*gradN_2');
% %         kk_13 = fac_Gauss'*diag(gradN_1*gradN_3');
% %         kk_14 = fac_Gauss'*diag(gradN_1*gradN_4');
% %         kk_15 = fac_Gauss'*diag(gradN_1*gradN_5');
% %         kk_16 = fac_Gauss'*diag(gradN_1*gradN_6');

% %         kk_23 = fac_Gauss'*diag(gradN_2*gradN_3');
% %         kk_24 = fac_Gauss'*diag(gradN_2*gradN_4');
% %         kk_25 = fac_Gauss'*diag(gradN_2*gradN_5');
% %         kk_26 = fac_Gauss'*diag(gradN_2*gradN_6');

% %         kk_34 = fac_Gauss'*diag(gradN_3*gradN_4');
% %         kk_35 = fac_Gauss'*diag(gradN_3*gradN_5');
% %         kk_36 = fac_Gauss'*diag(gradN_3*gradN_6');

% %         kk_45 = fac_Gauss'*diag(gradN_4*gradN_5');
% %         kk_46 = fac_Gauss'*diag(gradN_4*gradN_6');
% %         
% %         kk_56 = fac_Gauss'*diag(gradN_5*gradN_6');

        
        kk_loc = [kk_11 kk_12 kk_13 kk_14 kk_15 kk_16; ...
                  kk_12 kk_22 kk_23 kk_24 kk_25 kk_26; ...
                  kk_13 kk_23 kk_33 kk_34 kk_35 kk_36; ...
                  kk_14 kk_24 kk_34 kk_44 kk_45 kk_46; ...
                  kk_15 kk_25 kk_35 kk_45 kk_55 kk_56; ...
                  kk_16 kk_26 kk_36 kk_46 kk_56 kk_66];
                              
        % Sparse indexing
        num_el = numel(kk_loc);
        idx = num_el*(ii-1)+1:num_el*ii;    % block of 9 indices belonging to i-th triangle
        irow(idx) = repmat(tri_ii',3*N_order,1);    % row pointers
        jcol(idx) = reshape(repmat(tri_ii,3*N_order,1),num_el,1);    % column pointers
        
        KK_vals(idx) = kk_loc(:);    % nonzeros
        
    end
    
end

















% % coeffs = shape_f_norm(ii,:);
% % coeffs = reshape(coeffs,3*N_order,3*N_order)';
% % 
% % aa = coeffs(:,1);
% % bb = coeffs(:,2);
% % cc = coeffs(:,3);
% % dd = coeffs(:,4);
% % ee = coeffs(:,5);
% % ai = aa(1);
% % bi = bb(1);
% % ci = cc(1);
% % di = dd(1);
% % ei = ee(1);
% % aj = aa(1);
% % bj = bb(1);
% % cj = cc(1);
% % dj = dd(1);
% % ej = ee(1);
% % int_rz = @(r,z) ((2*ai*r+ci*z+di).*(2*aj*r+cj*z+dj) + (2*bi*z+ci*r+ei).*(2*bj*z+cj*r+ej))./(r-2)
% % 
% % P1 = [0 0];
% % P2 = [1 0];
% % P3 = [0 1];
% % 
% % xmin = min(nodes_Pla_ii(:,1));
% % xmax = max(nodes_Pla_ii(:,1));
% % ymin = min(nodes_Pla_ii(:,2));
% % ymax = max(nodes_Pla_ii(:,2));
% % I_rect = integral2(int_rz,xmin,xmax,ymin,ymax,'Method','iterated','AbsTol',0,'RelTol',1e-12);
% % 
% % m12 = (P2(2)-P1(2))/(P2(1)-P1(1));
% % xmin = min(nodes_Pla_ii([1 2],1));
% % xmax = max(nodes_Pla_ii([1 2],1));
% % ymin = min(nodes_Pla_ii([1 2],2));
% % ymax = @(x) (x-P1(1))*m12 + P1(2);
% % I_1 = integral2(int_rz,xmin,xmax,ymin,ymax,'Method','iterated','AbsTol',0,'RelTol',1e-12);
% % 
% % m23 = (P3(2)-P2(2))/(P3(1)-P2(1));
% % xmin = min(nodes_Pla_ii([2 3],1));
% % xmax = max(nodes_Pla_ii([2 3],1));
% % ymin = min(nodes_Pla_ii([2 3],2));
% % ymax = @(x) (x-P2(1))*m23 + P2(2);
% % I_2 = integral2(int_rz,xmin,xmax,ymin,ymax,'Method','iterated','AbsTol',0,'RelTol',1e-12);
% % 
% % m31 = (P1(2)-P3(2))/(P1(1)-P3(1));
% % xmin = min(nodes_Pla_ii([1 3],1));
% % xmax = max(nodes_Pla_ii([1 3],1));
% % ymin = min(nodes_Pla_ii([1 3],2));
% % ymax = @(x) (x-P3(1))*m31 + P3(2);
% % I_3 = integral2(int_rz,xmin,xmax,ymin,ymax,'Method','iterated','AbsTol',0,'RelTol',1e-12);
% % 
% % res_ref = I_rect - I_1 - I_2 - I_3
% % 
% % xmin = 0;
% % xmax = 1;
% % ymin = 0;
% % ymax = @(x) 1-x;
% % res_ref = integral2(int_rz,xmin,xmax,ymin,ymax,'Method','iterated','AbsTol',0,'RelTol',eps);
% % 
% % clear res
% % for kk = 2:49
% %     
% %     [ww,nodes_G,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,kk);
% %     
% %     rr_G_ii = nodes_G(:,1);
% %     zz_G_ii = nodes_G(:,2);
% %     
% %     gradN_1 = [2*aa(1)*rr_G_ii+cc(1)*zz_G_ii+dd(1)  2*bb(1)*zz_G_ii+cc(1)*rr_G_ii+ee(1)];
% %     
% %     fac_Gauss = ww.*(1./(rr_G_ii-2));
% %     
% %     % %             res(kk,:) = [kk n_Gauss fac_Gauss'*diag(gradN_1*gradN_1')];
% %     res = fac_Gauss'*diag(gradN_1*gradN_1');
% %     fprintf('%2i %3i %16.16f %e\n', kk, n_Gauss, res, abs(res_ref - res));
% %     
% % end









