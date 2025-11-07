function  [Source,solk] = fun_calcSourceFEM_adaptive(meshData_loc,solk,plaparameter,SETTINGS)



%%
mu0=4*pi*1.e-7;

ind_n_ONpla = meshData_loc.ind_n_ONpla;
ind_B = meshData_loc.ind_B;
ind_D = meshData_loc.ind_D;
ind_n_bc = meshData_loc.ind_n_bc;

KK_nobc = meshData_loc.KK_nobc;

tri             = meshData_loc.t;
N_order         = meshData_loc.shape_functions_N_order;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
n_Gauss         = meshData_loc.n_Gauss;
shape_functions = meshData_loc.shape_functions;
ww_Gauss = meshData_loc.ww_Gauss_pla;
nn = meshData_loc.nn;

rr_Gauss = P_Gauss(:,1);

ind_t_on_curve = meshData_loc.ind_t_on_curve;


%%

if SETTINGS.RUN_MEX_ROUTINE == false
    [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,solk.Psi,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,solk.Psi,N_order,P_Gauss,n_Gauss,shape_functions);
end

% %     [Psi_Gauss] = fun_evaluateFluxGaussPoints(meshData_loc,solk_fix.Psi);

% % Psi_k=Psi_Gauss;
% % Psi_axis=solk.Psi_axis;
% % Psi_Bpla=sum(solk.Psi_B)./numel(solk.Psi_B);
% % psi_bar=(Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux

Psi_k = Psi_Gauss;
Psi_axis = solk.Psi_axis;
Psi_Bpla = solk.Psi_B;
psi_bar = (Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux

solk.psi_bar=psi_bar;
solk.Psi_k=Psi_k;




%% Integration: triangles fully inside the plasma domain
ind_t_selected = setdiff(meshData_loc.ind_t_INpla,meshData_loc.ind_t_on_curve);

if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    Centroid = plaparameter.Centroid;
    Ipla =     plaparameter.Ipla;
    psibar =   plaparameter.psibar;
    FdF =      plaparameter.FdF;
    dP =      plaparameter.dP;
    
    psi_bar(psi_bar>1) = NaN;
    psi_bar(psi_bar<0) = NaN;
    
    fdfn=interp1(psibar,FdF,psi_bar,'linear');
    dpn=interp1(psibar,dP,psi_bar,'linear');
    gg=(rr_Gauss.*dpn+fdfn./(mu0*rr_Gauss)); %computation of nodes plasma current density
    gg=gg/2/pi;
    
    gg(isnan(gg)) = 0;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source_1] = fun_calcSourceFEM_v2(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        [Source_1] = fun_calcSourceFEM_v2_mex(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
        
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    
    R_0 =      plaparameter.R_0;
    beta_0 =   plaparameter.beta_0;
    alpha_M =  plaparameter.alpha_M;
    alpha_N =  plaparameter.alpha_N;
    Ipla =     plaparameter.Ipla;
    
    tmp_gg = ((1-psi_bar.^alpha_M).^alpha_N);
    tmp_gg(imag(tmp_gg) ~= 0 ) = abs(tmp_gg(imag(tmp_gg) ~= 0)); % occhio!!!
    aa_r = (rr_Gauss*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss);
    gg = tmp_gg.*aa_r;

    solk.aa_r_tot = aa_r;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source_1] = fun_calcSourceFEM_v2(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        [Source_1] = fun_calcSourceFEM_v2_mex(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
        
end

solk.gg_tot = gg;


%% Integration: triangles on the plasma domain

ind_t_selected = meshData_loc.ind_t_on_curve;

ind_n_vertices_plasma = meshData_loc.ind_n_vertices_plasma;
ind_n_plasma = meshData_loc.ind_n_INpla;
tri = meshData_loc.t;

omega_plasma = solk.Separatrix;

Psi_nodes = solk.Psi;

ind_t_Xp = solk.ind_t_Xp;
if isempty(ind_t_Xp)
    ind_t_Xp = NaN;
% %     ind_t_selected = ind_t_selected(~isnan(ind_t_selected));
end

tri_sel = tri(ind_t_selected,:);
nt = size(tri_sel,1);


nodes = meshData_loc.n;

% % triplot(tri_sel(:,1:3),nodes(:,1),nodes(:,2),'g')


Source_temp = zeros(nt,3*N_order);
Ind_temp = zeros(nt,3*N_order);

P1 = [0 0];
P2 = [1 0];
P3 = [0 1];

GAUSS_QUAD_DEGREE_PLA = SETTINGS.GAUSS_QUAD_DEGREE_PLA;

[ww_G_norm,P_G_norm,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,GAUSS_QUAD_DEGREE_PLA);


for jj = 1:nt
    
    ind_t_ii = ind_t_selected(jj);

    tri_ii = tri_sel(jj,:);
        
    % %         if 0
    % % % %             figure(ind_figure);
    % %             triplot(tri_ii(1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'k', 'linewidth' ,2); hold on
    % %             plot(omega_plasma(:,1),omega_plasma(:,2),'k')
    % %         end
        
        nodes_jj = nodes(tri_ii,:);
        Psi_nodes_jj = Psi_nodes(tri_ii);
        
        ind_n_vertices_inside_jj = ismember(tri_ii(1:3),ind_n_vertices_plasma);
        % %         ind_n_inside_jj = ismember(tri_ii,ind_n_plasma);
        
        L1 = nodes_jj([1:3 1],:)';
        L2 = omega_plasma([1:end 1],:)';
        P2 = InterX(L1,L2); P2 = P2';
        
        if size(P2,1) == 2 &&  ind_t_ii ~= ind_t_Xp % only one node inside
            type_of_case = 1; %one/two nodes inside, only one intersection
        elseif size(P2,1) > 2 && ind_t_ii ~= ind_t_Xp % more intersections, close to xpoint
            type_of_case = 2;
        elseif ind_t_ii == ind_t_Xp
            type_of_case = 3;
        end
        
        switch type_of_case
            case 1
                
                point_int_jj = nodes_jj(ind_n_vertices_inside_jj,:);
                
                % points on the boundary (vertices)
                L1 = nodes_jj([1:3 1],:)';
                L2 = omega_plasma([1:end 1],:)';
                P_intersection = InterX(L1,L2); P_intersection = P_intersection';
                
                poligono = [point_int_jj; P_intersection];
                poligono = fun_ordinapunti(poligono);
                
                centro = sum(poligono)/size(poligono,1);
                
                shape_f_jj = shape_functions(ind_t_ii,:);
                I_tri = 0;
                for kk=1:size(poligono,1)
                    if kk == size(poligono,1)
                        ind_sel = [kk 1];
                    else
                        ind_sel = [kk kk+1];
                    end
                    tri_kk = [poligono(ind_sel,:); centro];
                    tri_kk = fun_ordinapunti(tri_kk);
                    
                    P1 = tri_kk(1,:);
                    P2 = tri_kk(2,:);
                    P3 = tri_kk(3,:);
                    
% %                     plot(tri_kk([1:3 1],1),tri_kk([1:3 1],2),'-r', 'LineWidth',2)
% %                     plot(PP_G(:,1),PP_G(:,2),'ok', 'LineWidth',2)

                    [PP_G,ww_G] = fun_map_Gauss_point(P_G_norm,ww_G_norm,n_Gauss,P1,P2,P3);
                    rr_Gauss_kk = PP_G(:,1);
                    coeff = reshape(shape_f_jj,3*N_order,3*N_order)';
                    tmp_geo = [PP_G.^2 PP_G(:,1).*PP_G(:,2) ...
                        PP_G ones(n_Gauss,1)]';
                    
                    Psi_k_kk = (coeff*tmp_geo)'*Psi_nodes_jj;
                    psi_bar_jj = (Psi_k_kk-Psi_axis)/(Psi_Bpla-Psi_axis);
                    
                    if SETTINGS.J_PARAMETRIZATION_TYPE == 1
                        
                        psi_bar_jj(psi_bar>1) = NaN;
                        psi_bar_jj(psi_bar<0) = NaN;
                        
                        fdfn=interp1(psibar,FdF,psi_bar_jj,'linear');
                        dpn=interp1(psibar,dP,psi_bar_jj,'linear');
                        gg_kk=(rr_Gauss_kk.*dpn+fdfn./(mu0*rr_Gauss_kk)); %computation of nodes plasma current density
                        gg_kk=gg_kk/2/pi;
                        
                        gg_kk(isnan(gg_kk)) = 0;
                        
                    elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
                        
                        tmp_gg_kk = ((1-psi_bar_jj.^alpha_M).^alpha_N);
                        tmp_gg_kk(imag(tmp_gg_kk) ~= 0 ) = abs(tmp_gg_kk(imag(tmp_gg_kk) ~= 0)); % occhio!!!
                        aa_r_kk = (rr_Gauss_kk*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss_kk);
                        gg_kk = tmp_gg_kk.*aa_r_kk;
                        
                    end
                    
                    I_tri_kk = coeff*tmp_geo*(ww_G.*(gg_kk));

                    I_tri = I_tri + I_tri_kk;
                end
                
            case 2 % triangle close to the Xpoint (it sometimes exists, sometimes not)
                
% %                 if PLOT_FIGURES == 'T'
% %                     figure(ind_figure); triplot(tri_ii(1:3),nodes(:,1),nodes(:,2),'g', 'linewidth' ,2);
% %                 end
                
                point_int_jj = nodes_jj(ind_n_vertices_inside_jj,:);
                
                % points on the boundary (vertices)
                L1 = nodes_jj([1:3 1],:)';
                L2 = omega_plasma([1:end 1],:)';
                P_intersection = InterX(L1,L2); P_intersection = P_intersection';
                
                poligono = [point_int_jj; P_intersection];
                poligono = fun_ordinapunti(poligono);
                
                centro = sum(poligono)/size(poligono,1);
                
                shape_f_jj = shape_functions(ind_t_ii,:);
                
                I_tri = 0;
                for kk=1:size(poligono,1)
                    if kk == size(poligono,1)
                        ind_sel = [kk 1];
                    else
                        ind_sel = [kk kk+1];
                    end
                    tri_kk = [poligono(ind_sel,:); centro];
                    tri_kk = fun_ordinapunti(tri_kk);
                    
                    P1 = tri_kk(1,:);
                    P2 = tri_kk(2,:);
                    P3 = tri_kk(3,:);
                    
% %                     plot(tri_kk([1:3 1],1),tri_kk([1:3 1],2),'-r', 'LineWidth',2)
% %                     plot(PP_G(:,1),PP_G(:,2),'ok', 'LineWidth',2)

                    [PP_G,ww_G] = fun_map_Gauss_point(P_G_norm,ww_G_norm,n_Gauss,P1,P2,P3);
                    rr_Gauss_kk = PP_G(:,1);
                    coeff = reshape(shape_f_jj,3*N_order,3*N_order)';
                    tmp_geo = [PP_G.^2 PP_G(:,1).*PP_G(:,2) ...
                        PP_G ones(n_Gauss,1)]';
                    
                    Psi_k_kk = (coeff*tmp_geo)'*Psi_nodes_jj;
                    psi_bar_jj = (Psi_k_kk-Psi_axis)/(Psi_Bpla-Psi_axis);
                    
                    if SETTINGS.J_PARAMETRIZATION_TYPE == 1
                        
                        psi_bar_jj(psi_bar>1) = NaN;
                        psi_bar_jj(psi_bar<0) = NaN;
                        
                        fdfn=interp1(psibar,FdF,psi_bar_jj,'linear');
                        dpn=interp1(psibar,dP,psi_bar_jj,'linear');
                        gg_kk=(rr_Gauss_kk.*dpn+fdfn./(mu0*rr_Gauss_kk)); %computation of nodes plasma current density
                        gg_kk=gg_kk/2/pi;
                        
                        gg_kk(isnan(gg_kk)) = 0;
                        
                    elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
                        
                        tmp_gg_kk = ((1-psi_bar_jj.^alpha_M).^alpha_N);
                        tmp_gg_kk(imag(tmp_gg_kk) ~= 0 ) = abs(tmp_gg_kk(imag(tmp_gg_kk) ~= 0)); % occhio!!!
                        aa_r_kk = (rr_Gauss_kk*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss_kk);
                        gg_kk = tmp_gg_kk.*aa_r_kk;
                        
                    end
                    
                    I_tri_kk = coeff*tmp_geo*(ww_G.*(gg_kk));

                    I_tri = I_tri + I_tri_kk;
                end
                  
            case 3 % triangle with Xpoint
                
% %                 if PLOT_FIGURES == 'T'
% %                     figure(ind_figure); triplot(tri_ii(1:3),nodes(:,1),nodes(:,2),'g', 'linewidth' ,2);
% %                 end
                
                point_int_jj = nodes_jj(ind_n_vertices_inside_jj,:);
                
                % points on the boundary (vertices)
                L1 = nodes_jj([1:3 1],:)';
                L2 = omega_plasma([1:end 1],:)';
                P_intersection = InterX(L1,L2); P_intersection = P_intersection';
                
                poligono = [point_int_jj; P_intersection; Xp];
                poligono = fun_ordinapunti(poligono);
                
                centro = sum(poligono)/size(poligono,1);
                
                shape_f_jj = shape_functions(ind_t_ii,:);
                I_tri = 0;
                for kk=1:size(poligono,1)
                    if kk == size(poligono,1)
                        ind_sel = [kk 1];
                    else
                        ind_sel = [kk kk+1];
                    end
                    tri_kk = [poligono(ind_sel,:); centro];
                    tri_kk = fun_ordinapunti(tri_kk);
                    
                    P1 = tri_kk(1,:);
                    P2 = tri_kk(2,:);
                    P3 = tri_kk(3,:);
                    
% %                     plot(tri_kk([1:3 1],1),tri_kk([1:3 1],2),'-r', 'LineWidth',2)
% %                     plot(PP_G(:,1),PP_G(:,2),'ok', 'LineWidth',2)

                    [PP_G,ww_G] = fun_map_Gauss_point(P_G_norm,ww_G_norm,n_Gauss,P1,P2,P3);
                    rr_Gauss_kk = PP_G(:,1);
                    coeff = reshape(shape_f_jj,3*N_order,3*N_order)';
                    tmp_geo = [PP_G.^2 PP_G(:,1).*PP_G(:,2) ...
                        PP_G ones(n_Gauss,1)]';
                    
                    Psi_k_kk = (coeff*tmp_geo)'*Psi_nodes_jj;
                    psi_bar_jj = (Psi_k_kk-Psi_axis)/(Psi_Bpla-Psi_axis);
                    
                    if SETTINGS.J_PARAMETRIZATION_TYPE == 1
                        
                        psi_bar_jj(psi_bar>1) = NaN;
                        psi_bar_jj(psi_bar<0) = NaN;
                        
                        fdfn=interp1(psibar,FdF,psi_bar_jj,'linear');
                        dpn=interp1(psibar,dP,psi_bar_jj,'linear');
                        gg_kk=(rr_Gauss_kk.*dpn+fdfn./(mu0*rr_Gauss_kk)); %computation of nodes plasma current density
                        gg_kk=gg_kk/2/pi;
                        
                        gg_kk(isnan(gg_kk)) = 0;
                        
                    elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
                        
                        tmp_gg_kk = ((1-psi_bar_jj.^alpha_M).^alpha_N);
                        tmp_gg_kk(imag(tmp_gg_kk) ~= 0 ) = abs(tmp_gg_kk(imag(tmp_gg_kk) ~= 0)); % occhio!!!
                        aa_r_kk = (rr_Gauss_kk*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss_kk);
                        gg_kk = tmp_gg_kk.*aa_r_kk;
                        
                    end
                    
                    I_tri_kk = coeff*tmp_geo*(ww_G.*(gg_kk));

                    I_tri = I_tri + I_tri_kk;
                end
                
                
        end
        
        Ind_temp(jj,:) = tri_ii;
        Source_temp(jj,:) =  I_tri;
        
end
    

%
ind_vec = Ind_temp(:);
Source_vec = Source_temp(:);

Source_2 = accumarray(ind_vec,Source_vec,[nn 1]);


%%


Source = Source_1 + Source_2;

end

















