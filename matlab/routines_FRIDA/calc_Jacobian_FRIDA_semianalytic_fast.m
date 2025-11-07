%% R_psi

% DRpsi_Dpsi
dPsi = perturbation;

J_PARAMETRIZATION_TYPE = SETTINGS.J_PARAMETRIZATION_TYPE;

MatInd_nodes_tri = meshData_loc.MatInd_nodes_tri;
MatInd_nodes = meshData_loc.MatInd_nodes;
triangles = meshData_loc.t;
dim_ind_rows = 3*N_order*nnz(MatInd_nodes_tri(ind_n_INpla,:));
P_Gauss = meshData_loc.nodes_Gauss_pla;
n_Gauss = meshData_loc.n_Gauss;
Psi_nodes = solk.Psi;


if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    
    psibar = plaparameter.psibar;
    FdF = plaparameter.FdF;
    dP = plaparameter.dP;
    
% %     if numel(psibar) ~= 101
% %         
% %         psibar_new = linspace(psibar(1),psibar(end),101);
% %         FdF_new = interp1(psibar,FdF,psibar_new);
% %         dP_new = interp1(psibar,dP,psibar_new);
% %         
% % % %         figure
% % % %         subplot(1,2,1); plot(psibar,FdF,'k'); hold on;  plot(psibar_new,FdF_new,'or');
% % % %         subplot(1,2,2); plot(psibar,dP); hold on;  plot(psibar_new,dP_new,'or');
% %                 
% %     end
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        
        [ind_rows,ind_cols,dSource_Dpsi_vals] = fun_dSource_Dpsi_type_1(...
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
            mu0);
        
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        
        nn_size_max_par = SETTINGS.NN_SIZE_MAX_PARALLEL_JACOBIAN;
        [ind_rows,ind_cols,dSource_Dpsi_vals] = fun_dSource_Dpsi_type_1_par_mex(...
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
            mu0,...
            nn_size_max_par);
        
    end
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    
    alpha_M = plaparameter.alpha_M;
    alpha_N = plaparameter.alpha_N;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        
        [ind_rows,ind_cols,dSource_Dpsi_vals] = fun_dSource_Dpsi_type_2(...
            triangles,...
            N_order,...
            dim_ind_rows,...
            P_Gauss, ...
            MatInd_nodes_tri,...
            dPsi,...
            Psi_nodes,...
            Psi_Bpla,...
            Psi_axis,...
            alpha_M,...
            alpha_N,...
            aa_r_tot,...
            gg_tot,...
            ww_Gauss,...
            n_Gauss,...
            shape_functions,...
            ind_n_INpla);
        
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        
        nn_size_max_par = SETTINGS.NN_SIZE_MAX_PARALLEL_JACOBIAN;
        [ind_rows,ind_cols,dSource_Dpsi_vals] = fun_dSource_Dpsi_type_2_par_mex(...
            triangles,...
            N_order,...
            dim_ind_rows,...
            P_Gauss, ...
            MatInd_nodes_tri,...
            dPsi,...
            Psi_nodes,...
            Psi_Bpla,...
            Psi_axis,...
            alpha_M,...
            alpha_N,...
            aa_r_tot,...
            gg_tot,...
            ww_Gauss,...
            n_Gauss,...
            shape_functions,...
            ind_n_INpla,...
            nn_size_max_par);
        
    end
    
end

DSource_Dpsi = sparse(ind_rows,ind_cols,dSource_Dpsi_vals,nn,nn);
DSource_Dpsi_D = DSource_Dpsi(ind_D,ind_D);
Jac_source = -lambda_star_k*DSource_Dpsi_D;
DRpsi_Dpsi = KK_D + Jac_source;



%%% DRpsi_Dpsibc
DRpsi_Dpsibc = KK_bc;


%%% DRpsi_Dpsia
Dpsia = perturbation;
Psia_pert = Psi_axis + Dpsia;
psi_bar_pert=(Psi_k-Psia_pert)/(Psi_Bpla-Psia_pert); %normalized poloidal magnetic flux

if SETTINGS.RUN_MEX_ROUTINE == false
    [psi_bar_pert] = fun_evaluateFluxGaussPoints_v2(tri,psi_bar_pert,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [psi_bar_pert] = fun_evaluateFluxGaussPoints_v2_mex(tri,psi_bar_pert,N_order,P_Gauss,n_Gauss,shape_functions);
end


if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    
% %     psi_bar_pert(psi_bar_pert>1) = NaN;
    psi_bar_pert(psi_bar_pert<0) = 0;
    
    fdfn_bar_ii=interp1(psibar,FdF,psi_bar_pert,'linear');   %cubic interpolation
    dpn_bar_ii=interp1(psibar,dP,psi_bar_pert,'linear');     %cubic interpolation
    gg_tot_pert=(rr_Gauss.*dpn_bar_ii+fdfn_bar_ii./(mu0*rr_Gauss)); %computation of nodes plasma current density
    gg_tot_pert=gg_tot_pert/2/pi;
    ind_NAN = isnan(gg_tot_pert);
    gg_tot_pert(ind_NAN) = 0;
    Dgg_Dpsia = (gg_tot_pert - gg_tot)/Dpsia;
    Dgg_Dpsia(ind_NAN) = 0;
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    tmp_gg_pert = ((1-psi_bar_pert.^alpha_M).^alpha_N);
    tmp_gg_pert(imag(tmp_gg_pert) ~= 0) = abs(tmp_gg_pert(imag(tmp_gg_pert) ~= 0)); % occhio!!!
    gg_tot_pert = tmp_gg_pert.*aa_r_tot;
    Dgg_Dpsia = (gg_tot_pert - gg_tot)/Dpsia;
    
end

Dgg_Dpsia(meshData_loc.ind_G_out) = 0;

nn = meshData_loc.nn;
ind_t_selected = meshData_loc.ind_t_INpla;
ww_Gauss = meshData_loc.ww_Gauss_pla;

if SETTINGS.RUN_MEX_ROUTINE == false
    [DSource_Dpsia] = fun_calcSourceFEM_v2(tri,nn,Dgg_Dpsia,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
    [DSource_Dpsia] = fun_calcSourceFEM_v2_mex(tri,nn,Dgg_Dpsia,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
end

% % [DSource_Dpsia] = fun_calcSourceFEM(meshData_loc,index_loc,Dgg_Dpsia);
DSource_Dpsia_D = DSource_Dpsia(ind_D);
DRpsi_Dpsia = -lambda_star_k*DSource_Dpsia_D;


%%% DRpsi_Dpsib
Dpsib = perturbation;
Psi_Bpla_pert = Psi_Bpla + Dpsib;
psi_bar_pert=(Psi_k-Psi_axis)/(Psi_Bpla_pert-Psi_axis); %normalized poloidal magnetic flux

if SETTINGS.RUN_MEX_ROUTINE == false
    [psi_bar_pert] = fun_evaluateFluxGaussPoints_v2(tri,psi_bar_pert,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [psi_bar_pert] = fun_evaluateFluxGaussPoints_v2_mex(tri,psi_bar_pert,N_order,P_Gauss,n_Gauss,shape_functions);
end

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    
% %     psi_bar_pert(psi_bar_pert>1) = NaN;
    psi_bar_pert(psi_bar_pert<0) = 0;
    
    fdfn_bar_ii=interp1(psibar,FdF,psi_bar_pert,'linear');   %cubic interpolation
    dpn_bar_ii=interp1(psibar,dP,psi_bar_pert,'linear');     %cubic interpolation
    gg_tot_pert=(rr_Gauss.*dpn_bar_ii+fdfn_bar_ii./(mu0*rr_Gauss)); %computation of nodes plasma current density
    gg_tot_pert=gg_tot_pert/2/pi;
    ind_NAN = isnan(gg_tot_pert);
    gg_tot_pert(ind_NAN) = 0;
    Dgg_Dpsib = (gg_tot_pert - gg_tot)/Dpsib;
    Dgg_Dpsib(ind_NAN) = 0;
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    tmp_gg_pert = ((1-psi_bar_pert.^alpha_M).^alpha_N);
    tmp_gg_pert(imag(tmp_gg_pert) ~= 0) = abs(tmp_gg_pert(imag(tmp_gg_pert) ~= 0)); % occhio!!!
    gg_tot_pert = tmp_gg_pert.*aa_r_tot;
    Dgg_Dpsib = (gg_tot_pert - gg_tot)/Dpsib;
    
end

Dgg_Dpsib(meshData_loc.ind_G_out) = 0;


nn = meshData_loc.nn;
ind_t_selected = meshData_loc.ind_t_INpla;
ww_Gauss = meshData_loc.ww_Gauss_pla;

if SETTINGS.RUN_MEX_ROUTINE == false
    [DSource_Dpsib] = fun_calcSourceFEM_v2(tri,nn,Dgg_Dpsib,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
    [DSource_Dpsib] = fun_calcSourceFEM_v2_mex(tri,nn,Dgg_Dpsib,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
end

% % [DSource_Dpsib] = fun_calcSourceFEM(meshData_loc,index_loc,Dgg_Dpsib);
DSource_Dpsib_D = DSource_Dpsib(ind_D);
DRpsi_Dpsib = -lambda_star_k*DSource_Dpsib_D;


%%% DRpsi_Dlambda
DRpsi_Dlambda = -Source_D;


%% R_psibc

%%% DRpsibc_Dpsi
DRpsibc_Dpsi = KK_D_tilde;


%%% DRpsibc_Dpsibc
DRpsibc_Dpsibc = KK_bc_tilde;


%%% DRpsibc_Dpsia
DRpsibc_Dpsia = zeros(nn_B,1);


%%% DRpsibc_Dpsib
DRpsibc_Dpsib = zeros(nn_B,1);


%%% DRpsibc_Dlambda
DRpsibc_Dlambda = zeros(nn_B,1);


%% R_psia

%%% DRpsia_Dpsia
DRpsia_Dpsia = -1;

%%% DRpsia_Dpsib
DRpsia_Dpsib = 0;

%%% DRpsia_Dlambda
DRpsia_Dlambda = 0;


%% R_psib

%%% DRpsib_Dpsia
DRpsib_Dpsia = 0;

%%% DRpsib_Dpsib
DRpsib_Dpsib = -1;

%%% DRpsib_Dlambda
DRpsib_Dlambda = 0;


%% R_lambda

%%% DRlambda_Dpsi
DRlambda_Dpsi = (lambda_star_k*sum(DSource_Dpsi_D));
% %     DRlambda_Dpsi = lambda_k*dSource_Dpsi_D';


%%% DRlambda_Dpsibc
DRlambda_Dpsibc = zeros(1,nn_B);


%%% DRlambda_Dpsia
DRlambda_Dpsia = lambda_star_k*sum(DSource_Dpsia_D);


%%% DRlambda_Dpsib
DRlambda_Dpsib = lambda_star_k*sum(DSource_Dpsib_D);


%%% DRlambda_Dlambda
DRlambda_Dlambda = sum(Source_D);




