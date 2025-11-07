function [solk1] = fun_MagneticAxis(meshData_loc,solk1,SETTINGS)


Psi_Gauss = solk1.Psi_Gauss;
Grad_Psi = solk1.Grad_Gauss;

if meshData_loc.shape_functions_N_order == 1
    
    tri             = meshData_loc.t;
    Psi_nodes       = solk1.Psi;
    N_order         = meshData_loc.shape_functions_N_order;
    P_Gauss         = meshData_loc.nodes_Gauss_pla;
    n_Gauss         = meshData_loc.n_Gauss;
    shape_functions = meshData_loc.shape_functions;

    % %     if SETTINGS.RUN_MEX_ROUTINE == false
    % %         [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    % %     elseif SETTINGS.RUN_MEX_ROUTINE == true
    % %         [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    % %     end
    
    [Psi_axis,ind_max] = max(Psi_Gauss);

    rr_axis = P_Gauss(ind_max,1);
    zz_axis = P_Gauss(ind_max,2);

    ind_t_axis = ceil(ind_max/n_Gauss);

    tri_axis = tri(ind_t_axis,:);

    ind_n_axis = tri_axis;
    
    solk1.Psi_axis=Psi_axis;
    solk1.Axis_RR=rr_axis;
    solk1.Axis_ZZ=zz_axis;
    solk1.ind_n_axis=ind_n_axis;

    
elseif meshData_loc.shape_functions_N_order == 2
    
    %% Identify the triangle containing the axis)
    
    tri             = meshData_loc.t;
    Psi_nodes       = solk1.Psi;
    N_order         = meshData_loc.shape_functions_N_order;
    P_Gauss         = meshData_loc.nodes_Gauss_pla;
    n_Gauss         = meshData_loc.n_Gauss;
    shape_functions = meshData_loc.shape_functions;

    % %     if SETTINGS.RUN_MEX_ROUTINE == false
    % %         [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    % %     elseif SETTINGS.RUN_MEX_ROUTINE == true
    % %         [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    % %     end
    % %
    % %     if SETTINGS.RUN_MEX_ROUTINE == false
    % %         [Grad_Psi] = fun_calcGradientGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
    % %     elseif SETTINGS.RUN_MEX_ROUTINE == true
    % %         [Grad_Psi] = fun_calcGradientGauss_FEM_mex(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
    % %     end
    
% %     vec_ind_max = [0 0];
% %     vec_ind_max(1) = (ind_max);
% %     check_psi1 = nnz(solk1.Psi(tri(ceil(ind_max/n_Gauss),:)) <= Psi_axis);
% %     ind_max1 = ind_max;
% % 
% %     
% %     %%% check 
% %     [Psi_sort,ind_Psi_sort] = sort(Psi_Gauss,'descend');
% %     [norm_Grad_sort,~] = sort(fun_vecnorm(Grad_Psi));
% %         
% %     norm_Grad = fun_vecnorm(Grad_Psi(ind_Psi_sort,:));
% %     thresh = 15*norm_Grad_sort(1);
% %     check_psi2 = 0;
% %     
% %     if norm_Grad_star > thresh
% %         
% %         jj = 1;
% %         norm_grad_jj = norm_Grad(jj);
% %         while norm_grad_jj > thresh && jj<=numel(Psi_sort)
% %             jj = jj+1;
% %             norm_grad_jj = norm_Grad(jj);
% %         end
% %         
% %         Psi_axis = Psi_sort(jj);
% %         ind_max2 = ind_Psi_sort(jj);
% %         check_psi2 = nnz(solk1.Psi(tri(ceil(ind_max2/n_Gauss),:)) <= Psi_axis);
% %     end
% %     
% %     if check_psi1 >= check_psi2
% %         ind_max = ind_max1;
% %     else
% %         ind_max = ind_max2;        
% %     end
% %     
% %     % %     [Psi_axis,ind_max] = max(Psi_Gauss);
% %     Psi_axis = Psi_Gauss(ind_max);
% %     rr_axis = P_Gauss(ind_max,1);
% %     zz_axis = P_Gauss(ind_max,2);
% %     ind_t_axis = ceil(ind_max/n_Gauss);
% %     tri_axis = tri(ind_t_axis,:);
% %     ind_n_axis = tri_axis;


% %     figure
% %     rgrid=linspace(min(meshData_loc.n(:,1)),max(meshData_loc.n(:,1)),200);
% %     zgrid=linspace(min(meshData_loc.n(:,2)),max(meshData_loc.n(:,2)),200);
% %     [RR,ZZ] = meshgrid(rgrid,zgrid);
% %     contour(RR,ZZ,griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ),100);
% %     contour(RR,ZZ,griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ),[Psi_axis Psi_axis]);
% %     hold on;

    
    %% mew strategy (look for min grad and det(Hess)>0)
    
    [Hess_Psi] = fun_calcHessianGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
    
    det_Hess = Hess_Psi(:,1).*Hess_Psi(:,3) - Hess_Psi(:,2).^2;
    norm_Grad = fun_vecnorm(Grad_Psi);
    
    % point where to search
    ind_search = find([det_Hess>0 & Hess_Psi(:,1)<0]);



    [~,ind_sort] = sort(norm_Grad(ind_search));
    det_Hess_sort = det_Hess(ind_search(ind_sort));
    norm_Grad_sort = norm_Grad(ind_search(ind_sort));
    a_11_sort = Hess_Psi(ind_search(ind_sort),1);
    P_Gauss_sort = P_Gauss(ind_search(ind_sort),:);

    ind_max = ind_search(ind_sort(1));

    if isfield(solk1,'Separatrix_old')
        KEEP_LOOKING = true;
        ii = 0;
        while KEEP_LOOKING && ii < length(ind_sort)
            ii = ii + 1;
            in = inpolygon(P_Gauss_sort(ii,1),P_Gauss_sort(ii,2),solk1.Separatrix_old(:,1),solk1.Separatrix_old(:,2));
            % %             plot3(P_Gauss_sort(ii,1),P_Gauss_sort(ii,2), ...
            % %                 norm_Grad_sort(ii),'ob', 'LineWidth',2)
            % %             pause
    
            if in
                KEEP_LOOKING = false;
            end
        end
        ind_max = ind_search(ind_sort(ii));
    end

    Psi_axis = Psi_Gauss(ind_max);
    norm_Grad_axis = norm_Grad(ind_max);
    rr_axis = P_Gauss(ind_max,1);
    zz_axis = P_Gauss(ind_max,2);
    ind_t_axis = ceil(ind_max/n_Gauss);
    tri_axis = tri(ind_t_axis,:);
    ind_n_axis = tri_axis;
 
% %     figure
% %     plot(P_Gauss_sort(:,1),P_Gauss_sort(:,2),'*r', 'LineWidth',2); hold on;
% %     plot(P_Gauss(ind_max,1),P_Gauss(ind_max,2),'*b', 'LineWidth',2); hold on;
% %     plot(solk1.Separatrix_old(:,1),solk1.Separatrix_old(:,2),'k')
% %   
% % 
% %     figure
% %     plot3(P_Gauss_sort(:,1),P_Gauss_sort(:,2), ...
% %         norm_Grad_sort,'*r', 'LineWidth',2)
% %     hold on
% %     plot3(P_Gauss(ind_max,1),P_Gauss(ind_max,2), ...
% %         norm_Grad(ind_max),'co', 'LineWidth',2)


    %% mew strategy (look for min grad and det(Hess)>0)
    % %
    % %
    % %     if det_Hess_sort(1) > 0 && a_11_sort(1) < 0
    % %         ind_max = ind_sort(1);
    % %     else
    % %         % %         qq = double([det_Hess_sort>0 a_11_sort<0]);
    % %         % %         qqq = sum(qq,2);
    % %         % %         ii = min(find(qqq == 2));
    % %         % %         ind_max = ind_sort(ii);
    % %
    % %         ii = 1;
    % %         KEEP_LOOKING = true;
    % %         while KEEP_LOOKING && ii < length(det_Hess_sort)
    % %             ii = ii+1;
    % %             det_Hess_sort_ii = det_Hess_sort(ii);
    % %             a_11_sort_ii     = a_11_sort(ii);
    % %             if det_Hess_sort_ii > 0 && a_11_sort_ii < 0
    % %                 if isfield(solk1,'Separatrix')
    % %                     in = inpolygon(P_Gauss_sort(ii,1),P_Gauss_sort(ii,2),solk1.Separatrix(:,1),solk1.Separatrix(:,2));
    % %                     if in
    % %                         KEEP_LOOKING = false;
    % %                     end
    % %                 else
    % %                     KEEP_LOOKING = false;
    % %                 end
    % %             end
    % %             ind_max = ind_sort(ii);
    % %         end
    % %     end
    % %
    % %     Psi_axis = Psi_Gauss(ind_max);
    % %     norm_Grad_axis = norm_Grad(ind_max);
    % %     rr_axis = P_Gauss(ind_max,1);
    % %     zz_axis = P_Gauss(ind_max,2);
    % %     ind_t_axis = ceil(ind_max/n_Gauss);
    % %     tri_axis = tri(ind_t_axis,:);
    % %     ind_n_axis = tri_axis;
    % %
    % %
    % %         figure
    % %     rgrid=linspace(min(meshData_loc.n(:,1)),max(meshData_loc.n(:,1)),200);
    % %     zgrid=linspace(min(meshData_loc.n(:,2)),max(meshData_loc.n(:,2)),200);
    % %     [RR,ZZ] = meshgrid(rgrid,zgrid);
    % %         surf(RR,ZZ,griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ), 'EdgeColor','none');
    % %         hold on
    % %         plot3(P_Gauss(ind_max,1),P_Gauss(ind_max,2),Psi_Gauss(ind_max),'*r', 'LineWidth',2)
    % %
    % % figure
    % % % % contour(RR,ZZ,griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ), 100);
    % % hold on; axis equal
    % % plot(P_Gauss(ind_max,1),P_Gauss(ind_max,2),'*b', 'LineWidth',2)
    % % contour(RR,ZZ,griddata(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),norm_Grad,RR,ZZ), 100);
    % % % % contour(RR,ZZ,griddata(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),norm_Grad,RR,ZZ), [0 0]);
    % % plot(P_Gauss(det_Hess>0 & Hess_Psi(:,1) < 0,1),P_Gauss(det_Hess>0 & Hess_Psi(:,1) < 0,2),'*m', 'LineWidth',2)
    % % plot(P_Gauss(ind_max,1),P_Gauss(ind_max,2),'*b', 'LineWidth',2)
    % %
    % % figure
    % % plot3(P_Gauss(:,1),P_Gauss(:,2), ...
    % %     norm_Grad(:),'*k', 'LineWidth',2)
    % % hold on
    % % plot3(P_Gauss(det_Hess>0 & Hess_Psi(:,1) < 0,1),P_Gauss(det_Hess>0 & Hess_Psi(:,1) < 0,2), ...
    % %     norm_Grad(det_Hess>0 & Hess_Psi(:,1) < 0),'*r', 'LineWidth',2)
    % % plot3(P_Gauss(ind_max,1),P_Gauss(ind_max,2), ...
    % %     norm_Grad(ind_max),'co', 'LineWidth',2)

    %% Psi axis and weigths for delta_a (A(x)x=b)

    Psi_ii = Psi_nodes(tri_axis);
    
    P1_ii = meshData_loc.n(tri_axis(1),:);
    P2_ii = meshData_loc.n(tri_axis(2),:);
    P3_ii = meshData_loc.n(tri_axis(3),:);

% %     coeffs = meshData_loc.shape_functions_norm;
% %     coeffs = reshape(coeffs,3*N_order,3*N_order)';
% % 
% %     aa = coeffs(:,1);
% %     bb = coeffs(:,2);
% %     cc = coeffs(:,3);
% %     dd = coeffs(:,4);
% %     ee = coeffs(:,5);
% %     ff = coeffs(:,6);
% %     
% %     % mapping from arbitrary triangle [P1,P2,P3] to unit triangle [(0,0),(1,0),(0,1)]  
% %     T = [P2_ii' P3_ii'] - [P1_ii' P1_ii'];
% %     T_inv = [T(2,2) -T(1,2); -T(2,1) T(1,1)]/(T(1,1)*T(2,2) - T(1,2)*T(2,1));
% % 
% % 
% %     tmp = T_inv*([rr_axis zz_axis]' - P1_ii');
% %     
% %     rr_axis_norm = tmp(1);
% %     zz_axis_norm = tmp(2);
% % 
% %     delta_a_weights = (aa*rr_axis_norm^2 + bb*zz_axis_norm^2 + ...
% %         cc*rr_axis_norm*zz_axis_norm + dd*rr_axis_norm + ee*zz_axis_norm + ff)';
% %     
% %     delta_a_weights*Psi_ii - Psi_axis
    
    
    coeffs = meshData_loc.shape_functions(ind_t_axis,:);
    coeffs = reshape(coeffs,3*N_order,3*N_order)';
    
    aa = coeffs(:,1);
    bb = coeffs(:,2);
    cc = coeffs(:,3);
    dd = coeffs(:,4);
    ee = coeffs(:,5);
    ff = coeffs(:,6);
    
    delta_a_weights = (aa*rr_axis^2 + bb*zz_axis^2 + ...
        cc*rr_axis*zz_axis + dd*rr_axis + ee*zz_axis + ff)';

    % %     abs(delta_a_weights*Psi_ii - Psi_axis)
    if abs(delta_a_weights*Psi_ii - Psi_axis) > 1e-6
        warning('weights for Magnetic Axis computed inaccurately'); 
    end

% %     Psi_axis = delta_a_weights*Psi_ii;
    
    solk1.Psi_axis=Psi_axis;
    solk1.Axis_RR=rr_axis;
    solk1.Axis_ZZ=zz_axis;
    solk1.ind_n_axis=ind_n_axis;
    solk1.delta_a_weights = delta_a_weights;
    
    
    
end







end

