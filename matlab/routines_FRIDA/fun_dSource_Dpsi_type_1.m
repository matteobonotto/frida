function [ind_rows,ind_cols,dSource_Dpsi_vals] = fun_dSource_Dpsi_type_1(...
        triangles,...
        N_order,...
        dim_ind_rows,...
        P_Gauss, ...
        MatInd_nodes_tri,...
        dPsi,...
        Psi_nodes,...
        Psi_Bpla,...
        Psi_axis,...
        psibar,...
        FdF,...
        dP,...
        gg_tot,...
        ww_Gauss,...
        n_Gauss,...
        shape_functions,...
        ind_n_INpla,...
        mu0) %#codegen
    
    

%%
ind_rows = zeros(dim_ind_rows,1);
ind_cols = ind_rows;
dSource_Dpsi_vals = ind_rows;

rr_Gauss = P_Gauss(:,1);

%%
ind_start = 1;

if N_order == 1
    
    for ii = 1:length(ind_n_INpla)
        ii_node = ind_n_INpla(ii);
        
        tmp = MatInd_nodes_tri(ii_node,:);
        nnz_tmp = nnz(tmp);
        
        ii_tri_near = tmp(1:nnz_tmp);        
        tri = triangles(ii_tri_near,:);
        
        f_rz_nodes = Psi_nodes(tri);
        f_rz_nodes(tri==ii_node) = f_rz_nodes(tri==ii_node)+dPsi;
        f_nodes = f_rz_nodes;
        
        vec_dSource_Dpsi = zeros(numel(tri),1);
        vec_nodes = zeros(numel(tri),1);
        
        
        for jj = 1:size(tri,1)
            
            tri_jj = tri(jj,:);
            
            f_nodes_jj = f_nodes(jj,:);
            
            ind_t_jj = ii_tri_near(jj);
            
            coeffs = shape_functions(ind_t_jj,:);
            
            coeffs_rep = reshape(coeffs,3*N_order,numel(coeffs)/(3*N_order))';
            
            ii_tri_Gauss = (ind_t_jj-1)*n_Gauss+1:ind_t_jj*n_Gauss;
            
            gg_tot_ii = gg_tot(ii_tri_Gauss);
            Points_ii = P_Gauss(ii_tri_Gauss,:);
            
            tmp = [Points_ii ones(n_Gauss,1)]';
                       
            f_Gauss_ii = zeros(n_Gauss,1);
            
            for kk = 1:3*N_order
                aa = repmat(coeffs_rep(kk,1),n_Gauss,1);
                bb = repmat(coeffs_rep(kk,2),n_Gauss,1);
                cc = repmat(coeffs_rep(kk,3),n_Gauss,1);
                
                temp_geo =  aa.*Points_ii(:,1) + ...
                    bb.*Points_ii(:,2) + ...
                    cc;
                
                f_Gauss_tmp = f_nodes_jj(kk).*temp_geo;
                
                f_Gauss_ii = f_Gauss_ii + f_Gauss_tmp;
                
            end
            
% %             f_Gauss_ii = coeffs_rep(1,:)*tmp.*f_nodes_jj(1) + ...
% %                 coeffs_rep(2,:)*tmp.*f_nodes_jj(2) + ...
% %                 coeffs_rep(3,:)*tmp.*f_nodes_jj(3);
            
            Psi_k_ii=f_Gauss_ii(:);
            psi_bar_ii=(Psi_k_ii-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux
            
            psi_bar_ii(psi_bar_ii>1) = NaN;
            psi_bar_ii(psi_bar_ii<0) = NaN;
            
            tmp_interp = [];
            tmp_interp = interp1(psibar',[FdF' dP'],psi_bar_ii,SETTINGS.J_PARAMETRIZATION.INTERP1_METHOD); 
            
            fdfn_bar_ii = tmp_interp(:,1);
            dpn_bar_ii  = tmp_interp(:,2);
            
            gg_pert_psi=(rr_Gauss(ii_tri_Gauss).*dpn_bar_ii + fdfn_bar_ii./(mu0*rr_Gauss(ii_tri_Gauss)));
            gg_pert_psi=gg_pert_psi/2/pi;
            ind_NAN = isnan(gg_pert_psi);
            gg_pert_psi(ind_NAN) = 0;
            
            Dgg_Dpsi = (gg_pert_psi - gg_tot_ii)/dPsi;
            Dgg_Dpsi(ind_NAN) = 0;


            ww = ww_Gauss(ii_tri_Gauss,:);
            
            NN_1 = coeffs_rep(1,:)*tmp;
            NN_2 = coeffs_rep(2,:)*tmp;
            NN_3 = coeffs_rep(3,:)*tmp;
            
            tmp_ii = ww.*Dgg_Dpsi;
            
            DI_tri_1 = tmp_ii'*NN_1';
            DI_tri_2 = tmp_ii'*NN_2';
            DI_tri_3 = tmp_ii'*NN_3';
            
            DI_tri = [DI_tri_1; DI_tri_2; DI_tri_3];
            vec_dSource_Dpsi((jj-1)*3*N_order+1:jj*3*N_order) = DI_tri;
            vec_nodes((jj-1)*3*N_order+1:jj*3*N_order) = tri_jj(:);
            
        end
        
        ind_end = ind_start+numel(vec_dSource_Dpsi)-1;
        ind_rows(ind_start:ind_end) = vec_nodes;
        ind_cols(ind_start:ind_end) = ii_node;
        dSource_Dpsi_vals(ind_start:ind_end) = vec_dSource_Dpsi;
        
        ind_start = ind_end+1;
        
    end
    
    
elseif N_order == 2
    
    for ii = 1:length(ind_n_INpla)
        ii_node = ind_n_INpla(ii);
        
        tmp = MatInd_nodes_tri(ii_node,:);
        nnz_tmp = nnz(tmp);
        
        ii_tri_near = tmp(1:nnz_tmp);
        
        tri = triangles(ii_tri_near,:);
        
        f_rz_nodes = Psi_nodes(tri);
        f_rz_nodes(tri==ii_node) = f_rz_nodes(tri==ii_node)+dPsi;
        f_nodes = f_rz_nodes;
        
        vec_dSource_Dpsi = zeros(numel(tri),1);
        vec_nodes = zeros(numel(tri),1);
        
        for jj = 1:size(tri,1)
            
            tri_jj = tri(jj,:);
            
            f_nodes_jj = f_nodes(jj,:);
            
            ind_t_jj = ii_tri_near(jj);
            
            coeffs = shape_functions(ind_t_jj,:);
            
            coeffs_rep = reshape(coeffs,3*N_order,numel(coeffs)/(3*N_order))';
            
            ii_tri_Gauss = (ind_t_jj-1)*n_Gauss+1:ind_t_jj*n_Gauss;
            
            gg_tot_ii = gg_tot(ii_tri_Gauss);
            Points_ii = P_Gauss(ii_tri_Gauss,:);
            
            tmp = [Points_ii.^2 Points_ii(:,1).*Points_ii(:,2) Points_ii ones(n_Gauss,1)]';
                       
            f_Gauss_ii = zeros(n_Gauss,1);
            
            for kk = 1:3*N_order
                aa = repmat(coeffs_rep(kk,1),n_Gauss,1);
                bb = repmat(coeffs_rep(kk,2),n_Gauss,1);
                cc = repmat(coeffs_rep(kk,3),n_Gauss,1);
                dd = repmat(coeffs_rep(kk,4),n_Gauss,1);
                ee = repmat(coeffs_rep(kk,5),n_Gauss,1);
                ff = repmat(coeffs_rep(kk,6),n_Gauss,1);
                
                temp_geo =  aa.*Points_ii(:,1).^2 + ...
                    bb.*Points_ii(:,2).^2 + ...
                    cc.*Points_ii(:,1).*Points_ii(:,2) + ...
                    dd.*Points_ii(:,1) + ...
                    ee.*Points_ii(:,2) + ...
                    ff;
                
                f_Gauss_tmp = f_nodes_jj(kk).*temp_geo;
                
                f_Gauss_ii = f_Gauss_ii + f_Gauss_tmp;
                
            end

% %             f_Gauss_ii2 = [coeffs_rep(1,:)*tmp.*f_nodes_jj(1) + ...
% %                 coeffs_rep(2,:)*tmp.*f_nodes_jj(2) + ...
% %                 coeffs_rep(3,:)*tmp.*f_nodes_jj(3) + ...
% %                 coeffs_rep(4,:)*tmp.*f_nodes_jj(4) + ...
% %                 coeffs_rep(5,:)*tmp.*f_nodes_jj(5) + ...
% %                 coeffs_rep(6,:)*tmp.*f_nodes_jj(6)]';

            Psi_k_ii=f_Gauss_ii(:);
            psi_bar_ii=(Psi_k_ii-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux
            
            psi_bar_ii(psi_bar_ii>1) = NaN;
            psi_bar_ii(psi_bar_ii<0) = NaN;
            
            psi_bar_ii(psi_bar_ii>1) = NaN;
            psi_bar_ii(psi_bar_ii<0) = NaN;            

            tmp_interp = [];
            tmp_interp = interp1(psibar',[FdF' dP'],psi_bar_ii,'linear'); 
            
            fdfn_bar_ii = tmp_interp(:,1);
            dpn_bar_ii  = tmp_interp(:,2);

            gg_pert_psi=(rr_Gauss(ii_tri_Gauss).*dpn_bar_ii + fdfn_bar_ii./(mu0*rr_Gauss(ii_tri_Gauss)));
            gg_pert_psi=gg_pert_psi/2/pi;
            ind_NAN = isnan(gg_pert_psi);
            gg_pert_psi(ind_NAN) = 0;
            
            Dgg_Dpsi = (gg_pert_psi - gg_tot_ii)/dPsi;
            Dgg_Dpsi(ind_NAN) = 0;
            
            ww = ww_Gauss(ii_tri_Gauss,:);
            
            NN_1 = coeffs_rep(1,:)*tmp;
            NN_2 = coeffs_rep(2,:)*tmp;
            NN_3 = coeffs_rep(3,:)*tmp;
            NN_4 = coeffs_rep(4,:)*tmp;
            NN_5 = coeffs_rep(5,:)*tmp;
            NN_6 = coeffs_rep(6,:)*tmp;
            
            tmp_ii = ww.*Dgg_Dpsi;
            
            DI_tri_1 = tmp_ii'*NN_1';
            DI_tri_2 = tmp_ii'*NN_2';
            DI_tri_3 = tmp_ii'*NN_3';
            DI_tri_4 = tmp_ii'*NN_4';
            DI_tri_5 = tmp_ii'*NN_5';
            DI_tri_6 = tmp_ii'*NN_6';
            
            DI_tri = [DI_tri_1; DI_tri_2; DI_tri_3; DI_tri_4; DI_tri_5; DI_tri_6];
            vec_dSource_Dpsi((jj-1)*3*N_order+1:jj*3*N_order) = DI_tri;
            vec_nodes((jj-1)*3*N_order+1:jj*3*N_order) = tri_jj(:);
            
        end
        
        ind_end = ind_start+numel(vec_dSource_Dpsi)-1;
        ind_rows(ind_start:ind_end) = vec_nodes;
        ind_cols(ind_start:ind_end) = ii_node;
        dSource_Dpsi_vals(ind_start:ind_end) = vec_dSource_Dpsi;
        
        ind_start = ind_end+1;
        
    end
end



