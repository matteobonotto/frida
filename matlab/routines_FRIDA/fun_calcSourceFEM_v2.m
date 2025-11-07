function  [Source] = fun_calcSourceFEM_v2(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions)

%%

tri_sel = tri(ind_t_selected,:);
nt = size(tri_sel,1);

Source_temp = zeros(nt,3*N_order);
Ind_temp = zeros(nt,3*N_order);

if N_order == 1
    
    for ii = 1:nt %parfor nella versione mexata
        
        ind_t_ii = ind_t_selected(ii);
        tri_ii = tri_sel(ii,:);

        ii_tri_Gauss = (ind_t_ii-1)*n_Gauss+1:ind_t_ii*n_Gauss;
        
        ww = ww_Gauss(ii_tri_Gauss);
        nodes = P_Gauss(ii_tri_Gauss,:);
        
        shape_functions_ii = shape_functions(ind_t_ii,:);
        
        cc = reshape(shape_functions_ii,3*N_order,3*N_order)';
        
        tmp = [nodes ones(n_Gauss,1)]';
        
        NN_1 = cc(1,:)*tmp;
        NN_2 = cc(2,:)*tmp;
        NN_3 = cc(3,:)*tmp;
        
        gg_Gauss_ii = gg(ii_tri_Gauss);
        
        tmp_ii = ww.*gg_Gauss_ii;
        
        I_tri_1 = tmp_ii'*NN_1';
        I_tri_2 = tmp_ii'*NN_2';
        I_tri_3 = tmp_ii'*NN_3';
        
        I_tri = [I_tri_1; I_tri_2; I_tri_3]';
        
        Ind_temp(ii,:) = tri_ii;
        Source_temp(ii,:) =  I_tri;
              
        
    end
    
elseif N_order == 2
    
    for ii = 1:nt %parfor nella versione mexata
        
        ind_t_ii = ind_t_selected(ii);
        tri_ii = tri_sel(ii,:);

        ii_tri_Gauss = (ind_t_ii-1)*n_Gauss+1:ind_t_ii*n_Gauss;
        
        ww = ww_Gauss(ii_tri_Gauss);
        nodes = P_Gauss(ii_tri_Gauss,:);
        
        shape_functions_ii = shape_functions(ind_t_ii,:);
        
        cc = reshape(shape_functions_ii,3*N_order,3*N_order)';
        
        tmp = [nodes.^2 nodes(:,1).*nodes(:,2) nodes ones(n_Gauss,1)]';
        
        NN_1 = cc(1,:)*tmp;
        NN_2 = cc(2,:)*tmp;
        NN_3 = cc(3,:)*tmp;
        NN_4 = cc(4,:)*tmp;
        NN_5 = cc(5,:)*tmp;
        NN_6 = cc(6,:)*tmp;
        
% %         figure
% %         plot3(nodes(:,1),nodes(:,2),NN_1.*gg_Gauss_ii','o')
        
% %         ww = ww([4:6 1:3 7:12]);
        
% %         ww.'*(NN_1.*gg_Gauss_ii')';
        
        gg_Gauss_ii = gg(ii_tri_Gauss);
        
        tmp_ii = ww.*gg_Gauss_ii;
        
        I_tri_1 = tmp_ii'*NN_1';
        I_tri_2 = tmp_ii'*NN_2';
        I_tri_3 = tmp_ii'*NN_3';
        I_tri_4 = tmp_ii'*NN_4';
        I_tri_5 = tmp_ii'*NN_5';
        I_tri_6 = tmp_ii'*NN_6';
        
        I_tri = [I_tri_1; I_tri_2; I_tri_3; I_tri_4; I_tri_5; I_tri_6]';
        
        Ind_temp(ii,:) = tri_ii;
        Source_temp(ii,:) =  I_tri;
                
    end
    
end

ind_vec = Ind_temp(:);
Source_vec = Source_temp(:);

Source = accumarray(ind_vec,Source_vec,[nn 1]);

% % if size(Source,1) < nn
% %    
% %     Source = [Source; zeros(nn-size(Source,1),1)];
% %     
% % end


% % toc








% % figure
% % plot3(meshData.n(:,1),meshData.n(:,2),Sources,'.')
