function [OUTPUT_FRIDA]=run_FRIDA_evol_CrankNicolson_fast(SETTINGS)

%%

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')

fprintf('\n\n\n')
fprintf('   F  FFFFFFFF     R  RRRRRRRR     I  IIIIIIII     D  DDDDDDDD     A  AAAAAAAA\n')
fprintf('   FF  FFFFFFFF    RR  RRRRRRRR    II  IIIIIIII    DD  DDDDDDDD    AA  AAAAAAAA\n')
fprintf('   FFF  FFFFFFF    RRR  RRRRRRR    III  IIIIIII    DDD  DDDDDDD    AAA  AAAAAAA\n')
fprintf('   FFF             RRR      RRR        I  I        DDD      DDD    AAA      AAA\n')
fprintf('   FFF             RRR      RRR        II          DDD      DDD    AAA      AAA\n')
fprintf('   FFFFFFF         RRRRRRR  RRR        IIII        DDD      DDD    AAAAAAA  AAA\n')
fprintf('   FFFFFFFF        RRRRRRRR RRR        IIII        DDD      DDD    AAAAAAAA AAA\n')
fprintf('   FFFFFFFFF       RRRRRRRRRRRR        IIII        DDD      DDD    AAAAAAAAAAAA\n')
fprintf('   FFF             RRR   RRR           IIII        DDD      DDD    AAA      AAA\n')
fprintf('   FFF             RRR    RRR      IIIIIIIIIII     DDDDDDDDDDDD    AAA      AAA\n')
fprintf('   FFF             RRR     RRR     IIIIIIIIIIII    DDDDDDDDDDDD    AAA      AAA\n')
fprintf('    FF              RR       RR     IIIIIIIIIII    DDDDDDDDDDD      AA       AA\n')
fprintf('\n\n\n')

fprintf('   FRee-boundary Integro-Differential Axisymmetric (solver)\n')
fprintf('   M. Bonotto, D. Abate\n\n')
fprintf('   2020/2021\n')
fprintf('   ver. 3.0\n')


fprintf('\n\n\n')
fprintf('\n\n\n')

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')



%% Load Data
load(['INPUT_FRIDA_equil_' SETTINGS.filename_save],'plaparameter', 'Conductors');



%%

% % meshData=          INPUT_FRIDA.meshData;
% % meshData_loc=      INPUT_FRIDA.meshData_loc;
% % plaparameter=      INPUT_FRIDA.plaparameter;
% % Conductors=        INPUT_FRIDA.Conductors;
% % index=             INPUT_FRIDA.index;
% % index_loc=         INPUT_FRIDA.index_loc;


mu0=4*pi*1.e-7;

%% Initialization of SETTINGS and INPUT

run_initialize_SETTINGS




%% Preprocessing
if SETTINGS.PREPROC == true
    fprintf('\n\n\n')
    fprintf('Running Preprocessing \n')
    
    time_start_preproc = tic;
    
    % load geometry
    load(['INPUT_FRIDA_geo_' SETTINGS.filename_geo],'meshData', 'meshData_loc')
    
    % run preprocessing
% %     run_FRIDA_preprocessing
    run_FRIDA_preprocessing_evol
    
    fprintf(' \n');
    fprintf('--> Done Pre-processing!!! \n');
    
    time_end_preproc = toc(time_start_preproc);
    
    fprintf('    TOTAL Eelapsed time is %5.1f seconds \n',  time_end_preproc)
    
else
    fprintf('Skipping Preprocessing \n')
    
    load(['INPUT_FRIDA_geo_preproc_' SETTINGS.filename_geo], ...
        'meshData', ...
        'meshData_loc')
    
    fprintf('\n')
    
end



%% time for code execution

fprintf('\n\n\n')
fprintf('Starting execution of FRIDA code\n')

time_start_code = tic;


%%

Rplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,1))...
    max(meshData_loc.n(meshData_loc.ind_n_bc,1))]+[-0.3 0.3];
Zplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,2))...
    max(meshData_loc.n(meshData_loc.ind_n_bc,2))]+[-0.3 0.3];

rr=meshData_loc.n(:,1);
zz=meshData_loc.n(:,2);

rgrid=linspace(Rplot(1),Rplot(2),200);
zgrid=linspace(Zplot(1),Zplot(2),200);

[RR,ZZ] = meshgrid(rgrid,zgrid);

ind_figure = [];
if SETTINGS.FIGURES == true
    ind_figure=ceil(rand*100);
    figure(ind_figure)
    hold on;    grid on;   axis equal; colormap jet; box on
    hold on, xlabel('r [m]'), ylabel('z [m]');
    edges=pdemesh(meshData.n',meshData.e_interface',[]);
    % %     plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k')
    set(edges,'color','k')
    set(edges,'linewidth',1)
    axis([Rplot Zplot])
end

%%

fprintf('Equilibrium parameters \n \n')

disp(plaparameter)

fprintf('\n \n')



%% Matrix for Implicit Euler (nonlinear ODE solution)

h_step = SETTINGS.time_new - SETTINGS.time_old;

R_pas = meshData.R_pas;
L_pas = meshData.L_pas;

M_act_pas = meshData.M_act_pas;
M_CS_pas = meshData.M_CS_pas;

ind_act = (1:size(M_act_pas,2))';
ind_CS = (ind_act(end)+1:size(M_act_pas,2)+size(M_CS_pas,2))';

KONNAX_passive = meshData.KONNAX_passive;

R_ode = KONNAX_passive*R_pas*KONNAX_passive';
L_ode = KONNAX_passive*L_pas*KONNAX_passive';

n_pas = size(R_ode,1);

EE = eye(n_pas);
L_ode_inv = L_ode\EE;

M_pas = [M_act_pas M_CS_pas];
M_ode = KONNAX_passive*M_pas;

A = -R_ode*L_ode_inv;
C = L_ode_inv;

T_Ipla_eq_D = meshData.T_Ipla_eq_D;

B_act = R_ode*L_ode_inv*M_ode(:,ind_act);
B_pla = R_ode*L_ode_inv*M_ode(:,ind_CS)*T_Ipla_eq_D;
D_act = -L_ode_inv*M_ode(:,ind_act);
D_pla = -L_ode_inv*M_ode(:,ind_CS)*T_Ipla_eq_D;

meshData.A = A;
meshData.C = C;
meshData.B_act = B_act;
meshData.B_pla = B_pla;
meshData.D_act = D_act;
meshData.D_pla = D_pla;

meshData.L_ode_inv = L_ode_inv;


% For Eddy current equation
KK_D  = meshData_loc.KK_nobc(meshData_loc.ind_D,meshData_loc.ind_D);
KK_bc = meshData_loc.KK_nobc(meshData_loc.ind_D,meshData_loc.ind_B);


A_phic = (EE - h_step*A/2);
A_phic_inv = A_phic\EE;

A_phic_old = (EE + h_step*A/2);

Q_tilde = h_step*.5*A_phic_inv*B_pla*KK_D/mu0;
Q_tilde_bc = h_step*.5*A_phic_inv*B_pla*KK_bc/mu0;


GG_tilde_bc = meshData.Green_pas_Bc*KONNAX_passive'/2/pi;

meshData.GG_tilde_bc = GG_tilde_bc;

P_tilde = -GG_tilde_bc*( C*Q_tilde + D_pla*KK_D/mu0);
P_tilde_bc = -GG_tilde_bc*( C*Q_tilde_bc + D_pla*KK_bc/mu0);


meshData.h_step = h_step;
meshData.A_phic = A_phic;
meshData.A_phic_inv = A_phic_inv;
meshData.A_phic_old = A_phic_old;
meshData.P_tilde = P_tilde;
meshData.P_tilde_bc = P_tilde_bc;
meshData.Q_tilde = Q_tilde;
meshData.Q_tilde_bc = Q_tilde_bc;




%%
% Vacuum solution
% % fprintf('Computing vacuum solution ... \n \n')
tic
% % [solk0] = fun_Vacuum_FRIDA(meshData,meshData_loc,Conductors,SETTINGS);
[solk0] = fun_Vacuum_FRIDA_evol(meshData,meshData_loc,Conductors,SETTINGS);
toc
fprintf('\n')

%
if SETTINGS.FIGURES_DEBUG == true
    figure
    subplot(1,2,1)
    edges=pdemesh(meshData.n',meshData.e_interface',[]);
    set(edges,'color','k')
    set(edges,'linewidth',1)
    axis equal;
    axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
    PSI = griddata(rr,zz,solk0.Psi_MS,RR,ZZ);
    contour(RR,ZZ,PSI,50);
% %     plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'.k')
    colormap('cool');
    title('Vacuum solution')
    colorbar vert
end

if SETTINGS.FIRST_SOLUTION == false    
    
    % First flux map
    fprintf('First flux map ... \n')
    tic
    [solk0] = fun_firstfluxmap(solk0,meshData,meshData_loc,plaparameter,SETTINGS);
    toc
    fprintf('\n')
    
    %
    if SETTINGS.FIGURES_DEBUG == true
        % %     figure
        subplot(1,2,2);
        edges=pdemesh(meshData.n',meshData.e_interface',[]);
        set(edges,'color','k')
        set(edges,'linewidth',1)
        axis equal;
        axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
        PSI = griddata(rr,zz,solk0.Psi,RR,ZZ);
        contour(RR,ZZ,PSI,30);
% %         plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k')
% %         colormap('cool');
        title('First solution')
        colorbar vert
        
    end
    
    % First Magnetic Axis
    fprintf('First Magnetic Axis ... \n')

    [solk0.Psi_axis,ind_axis]=max(solk0.Psi);
    
    solk0.ind_n_axis = ind_axis;
    solk0.delta_a_weights = 1;
    solk0.Axis_RR = solk0.r0;
    solk0.Axis_ZZ = solk0.z0;
    
    toc
    fprintf('\n')
        
elseif SETTINGS.FIRST_SOLUTION == true
    
    xx_k = SETTINGS.xx_0;
    
    ind_D = meshData_loc.ind_D; 
    ind_B = meshData_loc.ind_B; 
    
    nn_D = numel(ind_D);
    nn_B = numel(ind_B);
    
    temp_Psi          = xx_k(1:meshData_loc.nn);
    solk0.Psi         = zeros(meshData_loc.nn,1);
    solk0.Psi(ind_D)  = temp_Psi(1:nn_D);
    solk0.Psi(ind_B)  = temp_Psi(nn_D+1:nn_D+nn_B);
% %     solk0.Psi_axis   = xx_k(meshData_loc.nn+1);
% %     solk0.Psi_Bpla   = xx_k(meshData_loc.nn+2);
    solk0.lambda_star = xx_k(meshData_loc.nn+3);
    solk0.lambda      = solk0.lambda_star/mu0;
    
    % First Magnetic Axis
    fprintf('First Magnetic Axis ... \n')
    [solk0]=fun_MagneticAxis(meshData_loc,solk0,SETTINGS);
    
    
    % Compute fluxes at Gauss points
    fprintf('Compute fluxes at Gauss points ... \n')
    tic
    Psi_nodes        = solk0.Psi;
    tri              = meshData_loc.t;
    shape_functions  = meshData_loc.shape_functions;
    N_order          = meshData_loc.shape_functions_N_order;
    n_Gauss          = meshData_loc.n_Gauss;
    P_Gauss          = meshData_loc.nodes_Gauss_pla;

    if SETTINGS.RUN_MEX_ROUTINE == false
        [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    end
    solk0.Psi_Gauss = Psi_Gauss;
    toc; fprintf('\n');
    
        
    % First XPoint/Limiter and separatrix
    fprintf('First XPoint/Limiter and separatrix \n'); tic
    [solk0]=fun_Separatrix(meshData_loc,solk0,SETTINGS);
    toc; fprintf('\n');
    
    fprintf('Nodes inside plasma and on plasma boundary \n'); tic
    % Nodes inside plasma and on plasma boundary
    [meshData_loc]=fun_UpdateNodes_fast(meshData_loc,solk0);
    toc; fprintf('\n');
    
    
    % Complute plasma current
    fprintf('Complute plasma current \n'); tic
    
% %     if SETTINGS.J_PARAMETRIZATION_TYPE == 1
% %         Centroid = plaparameter.Centroid;
% %         Ipla =     plaparameter.Ipla;
% %         psibar =   plaparameter.psibar;
% %         FdF =      plaparameter.FdF;
% %         dP =      plaparameter.dP;
% %         
% %     elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
% %         R_0 =      plaparameter.R_0;
% %         beta_0 =   plaparameter.beta_0;
% %         alpha_M =  plaparameter.alpha_M;
% %         alpha_N =  plaparameter.alpha_N;
% %         Ipla =     plaparameter.Ipla;
% %     end
% %     
% %     rr_Gauss        = meshData_loc.nodes_Gauss_pla(:,1);
% % 
% %     [~,ind] = max(solk0.Psi(solk0.ind_n_axis));
% %     solk0.Psi_axis_node = solk0.Psi(solk0.ind_n_axis(ind));
% % 
% %     Psi_k    = solk0.Psi;
% %     Psi_axis = solk0.Psi_axis_node;
% %     Psi_Bpla = sum(solk0.Psi_Bpla)./numel(solk0.Psi_Bpla);
% %     psi_bar  = (Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux
% %     
% %     if SETTINGS.RUN_MEX_ROUTINE == false
% %         [psi_bar] = fun_evaluateFluxGaussPoints_v2(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
% %     elseif SETTINGS.RUN_MEX_ROUTINE == true
% %         [psi_bar] = fun_evaluateFluxGaussPoints_v2_mex(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
% %     end
% %     solk0.psi_bar = psi_bar;
% %     
% %     if SETTINGS.J_PARAMETRIZATION_TYPE == 1
% %         
% %         % %     psi_bar(psi_bar>1) = NaN;
% %         psi_bar(psi_bar<0) = 0;
% %         
% %         fdfn=interp1(psibar,FdF,psi_bar,'linear');
% %         dpn=interp1(psibar,dP,psi_bar,'linear');
% %         gg_tot=(rr_Gauss.*dpn+fdfn./(mu0*rr_Gauss)); %computation of nodes plasma current density
% %         gg_tot=gg_tot/2/pi;
% %         
% %         qq = isnan(gg_tot);
% %         gg_tot(isnan(gg_tot)) = 0;
% %         
% %         nn = meshData_loc.nn;
% %         ind_t_selected = meshData_loc.ind_t_INpla;
% %         ww_Gauss = meshData_loc.ww_Gauss_pla;
% %         
% %         if SETTINGS.RUN_MEX_ROUTINE == false
% %             [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
% %         elseif SETTINGS.RUN_MEX_ROUTINE == true
% %             if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
% %             [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
% %         end
% %         
% %         % %     figure;  plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),Source,'o')
% %         
% %     elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
% %         
% %         tmp_gg = ((1-psi_bar.^alpha_M).^alpha_N);
% %         tmp_gg(imag(tmp_gg) ~= 0) = abs(tmp_gg(imag(tmp_gg) ~= 0)); % occhio!!!
% %         aa_r_tot = (rr_Gauss*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss);
% %         gg_tot = tmp_gg.*aa_r_tot;
% %         
% %         nn = meshData_loc.nn;
% %         ind_t_selected = meshData_loc.ind_t_INpla;
% %         ww_Gauss = meshData_loc.ww_Gauss_pla;
% %         
% %         if SETTINGS.RUN_MEX_ROUTINE == false
% %             [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
% %         elseif SETTINGS.RUN_MEX_ROUTINE == true
% %             if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
% %             [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
% %         end
% %         
% %     end
% %     
% %    
% %     Iphi = Source*solk0.lambda;
% %     solk0.Iphi_old = Iphi;
    
    
    solk0.Iphi_old = (KK_D*SETTINGS.xx_0(1:nn_D) + KK_bc*SETTINGS.xx_0(nn_D+1:nn_D+nn_B))/mu0;
    
    toc; fprintf('\n');

    
    fprintf('\n')
    
end




% Passive currents 
solk0.ipass = SETTINGS.i_pas_old;
solk0.phi_pas_old = SETTINGS.phi_pas_old; 
solk0.i_pas_old = meshData.KONNAX_passive*SETTINGS.i_pas_old;
solk0.psi_bc_pas_old = SETTINGS.psi_bc_pas_old;


% Compute fluxes at Gauss points
fprintf('Compute fluxes at Gauss points ... \n')
tic
Psi_nodes        = solk0.Psi;
tri              = meshData_loc.t;
shape_functions  = meshData_loc.shape_functions;
N_order          = meshData_loc.shape_functions_N_order;
n_Gauss          = meshData_loc.n_Gauss;
P_Gauss          = meshData_loc.nodes_Gauss_pla;

if SETTINGS.RUN_MEX_ROUTINE == false
    [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
end
solk0.Psi_Gauss = Psi_Gauss;
toc


% First XPoint/Limiter and separatrix
fprintf('First XPoint/Limiter and separatrix ... \n')
tic
[solk0] = fun_Separatrix(meshData_loc,solk0,SETTINGS);

if SETTINGS.FIGURES == true
    % %     contour(RR,ZZ,PSI,[min(solk0.Psi_B) max(solk0.Psi_B)],'r','linewidth',1);
    plot(solk0.Separatrix(:,1),solk0.Separatrix(:,2),'k-')
    plot(solk0.Axis_RR,solk0.Axis_ZZ,'ok')
    % %     plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
    drawnow
end
toc

% Find nodes inside plasma and on plasma boundary
fprintf('Find nodes inside plasma and on plasma boundary ... \n')
tic
[meshData_loc]=fun_UpdateNodes_fast(meshData_loc,solk0);
% [index_loc]=fun_UpdateNodes(meshData_loc,solk0,index_loc);
% % plot(meshData_loc.n(index_loc.ind_n_ONpla,1),meshData_loc.n(index_loc.ind_n_ONpla,2),'*r')
% % plot(meshData_loc.n(index_loc.ind_n_INpla,1),meshData_loc.n(index_loc.ind_n_INpla,2),'ob')
toc
fprintf('\n')

pause(.0001)

% Smoothing of profiles

if SETTINGS.J_PARAMETRIZATION_TYPE == true
    
    if SETTINGS.SMOOTH_FDFPROFILE == true
        
        fprintf('Smoothing of FdF profile\n')
        
        npt_smoothing = SETTINGS.SMOOTH_FDFPROFILENPT;
        FdF_smooth = smoothdata(plaparameter.FdF,'sgolay',npt_smoothing);
        FdF_smooth(end) = 0;
        
        figure
        plot(plaparameter.psibar,plaparameter.FdF,'r-'); hold on;
        plot(plaparameter.psibar,FdF_smooth,'b--');
        title('FdF smoothing'); legend('Original', 'Smoothed')
        plaparameter.FdF = FdF_smooth;
        
    end
    
    if SETTINGS.SMOOTH_DPPROFILE == true
        
        fprintf('Smoothing of dP profile\n')
        
        npt_smoothing = SETTINGS.SMOOTH_DPPROFILENPT;
        dP_smooth = smoothdata(plaparameter.dP,'sgolay',npt_smoothing);
        dP_smooth(end) = 0;
        
        figure
        plot(plaparameter.psibar,plaparameter.dP,'r-'); hold on;
        plot(plaparameter.psibar,dP_smooth,'b--');
        title('dP smoothing'); legend('Original', 'Smoothed')
        plaparameter.dP = dP_smooth;
        
    end
    
end


if SETTINGS.FIRST_SOLUTION == false% %     Psi_k=solk0.Psi;
% %     solk0.Psi_k=Psi_k;

% %     xx_k=[solk0.Psi; solk0.Psi_axis; solk0.Psi_B;  0];

    Psi_0_D = solk0.Psi(meshData_loc.ind_D);
    Psi_0_B = solk0.Psi(meshData_loc.ind_B);
    xx_k = [Psi_0_D; ...
        Psi_0_B; ...
        solk0.Psi_axis; ...
        solk0.Psi_B; ...
        1]; 
    
end


solk0.ipass = SETTINGS.i_pas_old;
solk0.psi_pas_old = SETTINGS.psi_pas_old; 
solk0.i_pas_old = meshData.KONNAX_passive*SETTINGS.i_pas_old;
solk0.psi_bc_pas_old = SETTINGS.psi_bc_pas_old;


%% Iterative solution of FREE Boundary problem (FreeB)

SOLVER=SETTINGS.SOLVER;

nn=meshData_loc.nn;

% % solk0.lambda = xx_k(nn+3);
solk0.xx_k = xx_k;

% CICLO ITERATIVO
TOLL               = SETTINGS.TOLL;
SCARTO             = SETTINGS.SCARTO;
ii_freeB_max       = SETTINGS.ii_freeB_max;
RESIDUO_ii         = 100;
RESIDUO_norm_ii    = 100;
SCARTO_norm_ii     = 100;
ii_freeB           = 1;
scarto_tot         = [];
residuo_tot        = [];
residuo_norm_tot   = [];
time_per_iteration = [];
time_solver        = [];

%
solk0.residuo      = RESIDUO_ii;
solk0.residuo_norm = RESIDUO_norm_ii;
solk0.error_Ipla   = 0;
solk0.RR_psi_k_old = 0;

solk = solk0;



while RESIDUO_norm_ii > TOLL && SCARTO_norm_ii > SCARTO && ii_freeB <= ii_freeB_max
    
    time_start_iteration = tic;
    fprintf('\n');
    fprintf('Free Boundary iteration %i \n', ii_freeB);
    fprintf('\n');
    
    JP.faces           = meshData_loc.t(meshData_loc.ind_t_InFW,1:3);
    JP.vertices        = meshData_loc.n;
    JP.facevertexcdata = zeros(nn,1);
    
    %% Solver
    
    INPUT_solver.solk            = solk;
    INPUT_solver.meshData        = meshData;
    INPUT_solver.meshData_loc    = meshData_loc;
    INPUT_solver.plaparameter    = plaparameter;
    INPUT_solver.Conductors      = Conductors;
    INPUT_solver.SETTINGS        = SETTINGS;
    INPUT_solver.ii_freeB        = ii_freeB;
    INPUT_solver.RESIDUO_norm_ii = RESIDUO_norm_ii;
    
    
    
    switch SOLVER
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'PICARD'
% %             [solk1]=fun_Picard_FRIDA(INPUT_solver);
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'NR' % Newtok-Raphson
            [solk1]=fun_NR_FRIDA_evol_CrankNicolson_fast(INPUT_solver);
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'NK' % Newtok-Krylov
            [solk1]=fun_NK_FRIDA_evol_CrankNicolson(INPUT_solver);  %  Newtok-Krylov scaled lambda
            
    end
    
    
    RESIDUO_ii      = norm(solk1.h_k);
    RESIDUO_norm_ii = norm(solk1.h_k)/norm([solk1.xx_k(1:end-1); mu0*solk1.xx_k(end)]);
    SCARTO_norm_ii  = norm(solk1.xx_k-solk.xx_k)/norm(solk1.xx_k);
    
    solk1.residuo = RESIDUO_ii;
    solk1.residuo_norm = RESIDUO_norm_ii;
    
    fprintf('\n');
    
    %% Update quantities
    
    Psi_nodes        = solk1.Psi;
    tri              = meshData_loc.t;
    shape_functions  = meshData_loc.shape_functions;
    N_order          = meshData_loc.shape_functions_N_order;
    n_Gauss          = meshData_loc.n_Gauss;
    P_Gauss          = meshData_loc.nodes_Gauss_pla;

    
    % Compute Gradient of Psi (at Gauss points)
% %     fprintf('Compute Gradient of Psi (on nodes) ... \n')   
% %     tic
% %     if SETTINGS.RUN_MEX_ROUTINE == false
% %         [Grad_Gauss] = fun_calcGradientGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
% %     elseif SETTINGS.RUN_MEX_ROUTINE == true
% %         [Grad_Gauss] = fun_calcGradientGauss_FEM_mex(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
% %     end
% %     solk1.Grad_Gauss = Grad_Gauss;
% %     toc; fprintf('\n');

    
    % Compute fluxes at Gauss points
    fprintf('Compute fluxes at Gauss points ... \n')
    tic
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,Psi_nodes,N_order,P_Gauss,n_Gauss,shape_functions);
    end
    solk1.Psi_Gauss = Psi_Gauss;
    toc; fprintf('\n');


    % Update Magnetic Axis
    fprintf('Update Magnetic Axis \n'); tic
    [solk1] = fun_MagneticAxis(meshData_loc,solk1,SETTINGS);
    toc; fprintf('\n');
    
    % Update XPoint/Limiter and separatrix
    fprintf('Update XPoint/Limiter and separatrix \n'); tic
    [solk1] = fun_Separatrix(meshData_loc,solk1,SETTINGS);
    toc; fprintf('\n');
    
    % Update Nodes inside plasma and on plasma boundary
    fprintf('Update Nodes inside plasma and on plasma boundary \n'); tic
    [meshData_loc] = fun_UpdateNodes_fast(meshData_loc,solk1);
    toc; fprintf('\n');
    
    % Update Centroid
    fprintf('Update Centroid \n'); tic
    % %     [solk1] = fun_Centroid(meshData_loc,solk1);
    [r_p,z_p]=fun_Centroid_def(meshData_loc,solk1);
    solk1.Centroid_RR = r_p;
    solk1.Centroid_ZZ = z_p;
    toc; fprintf('\n');
    

    %% Plot results
    
    if SETTINGS.FIGURES_DEBUG == true
        figure
        subplot(1,3,1)
        edges=pdemesh(meshData.n',meshData.e_interface',[]);
        set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
        axis([Rplot Zplot]),  hold on,
        PSI = griddata(rr,zz,solk1.Psi,RR,ZZ);
        contour(RR,ZZ,PSI,100);
        colormap('cool'); colorbar vert;
        % %         plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
        plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-k');
        
        subplot(1,3,2)
        F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),solk1.Jphi);
        J_interp = zeros(meshData_loc.nn,1);
        J_interp(meshData_loc.ind_n_INpla) = F(meshData_loc.n(meshData_loc.ind_n_INpla,1),meshData_loc.n(meshData_loc.ind_n_INpla,2));
        JP.facevertexcdata=J_interp;
        hh=patch(JP,'facecolor','interp','edgecolor','none');
        axis equal;  hold on; colorbar vert; % title('j_\phi');
        axis([Rplot Zplot]),  hold on,
        axis square;  hold on;  colorbar vert; box on;
        edges=pdemesh(meshData.n',meshData.e_interface',[]);
        set(edges,'color','k'); set(edges,'linewidth',1)
        axis equal; axis([Rplot Zplot]),  hold on,
        % %         plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
        plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-k');
        
        subplot(1,3,3)
        H_K=JP;
        H_K.facevertexcdata=solk1.h_k(1:end-3);
        hh=patch(H_K,'facecolor','interp','edgecolor','none');
        axis square;  hold on;  colorbar vert; box on;
        edges=pdemesh(meshData.n',meshData.e_interface',[]);
        set(edges,'color','k'); set(edges,'linewidth',1)
        axis equal; axis([Rplot Zplot]),  hold on,
        % %         plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
        plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-k');
        
    end
    
    if SETTINGS.FIGURES == true
        h_psi = figure(ind_figure);
        [RR,ZZ] = meshgrid(rgrid,zgrid);
        rr = meshData_loc.n(:,1);
        zz = meshData_loc.n(:,2);
        PSI = griddata(rr,zz,solk1.Psi,RR,ZZ);
        % %         contour(RR,ZZ,PSI,100);
        % %         hold on
        % %         contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'k');
        plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-k');
        % %         plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
        plot(solk1.Axis_RR,solk1.Axis_ZZ,'ko')
        plot(solk1.XP_RR,solk1.XP_ZZ,'ko')
        drawnow
        
    end
    
    % %     if ii_freeB == 20
    % %         pause
    % %     end
    
    %%
    if SETTINGS.FIRST_SOLUTION == false    
        solk1.Psi_MS  = solk0.Psi_MS;
    end
    
    solk1.Psi_loc     = solk.Psi;
    solk1.phi_pas_old = solk.phi_pas_old;
    solk1.Iphi_old    = solk.Iphi_old;
% %     solk1.Iphi_old    = solk1.Iphi(ind_D);

    
    solk = solk1;
    
    fprintf(' \n')
    fprintf('RESIDUO                     %1.4e \n', RESIDUO_ii);
    fprintf('RESIDUO NORM                %1.4e \n', RESIDUO_norm_ii);
    fprintf('SCARTO                      %1.4e \n', SCARTO_norm_ii);
    fprintf('Check total plasma current  %10.10e \n', solk1.Ipla_check);
    
    time_end_iteration = toc(time_start_iteration);
    fprintf('Time per iteration is %5.6f seconds \n',  time_end_iteration)
    fprintf('\n \n');
    
    scarto_tot         = [scarto_tot; SCARTO_norm_ii];
    residuo_tot        = [residuo_tot; RESIDUO_ii];
    residuo_norm_tot   = [residuo_norm_tot; RESIDUO_norm_ii];
    time_per_iteration = [time_per_iteration; time_end_iteration];
    time_solver        = [time_solver; solk1.time_solver];
    
    ii_freeB = ii_freeB+1;
    
    
end

%% Impose zero current density outside separatrix (only for plotting purposes)
Boundary = solk1.Separatrix;
nodes_Gauss_pla = meshData_loc.nodes_Gauss_pla;
ind_Gauss_INpla=find(inpolygon(nodes_Gauss_pla(:,1),nodes_Gauss_pla(:,2),Boundary(:,1),Boundary(:,2)));

Jpla = zeros(size(nodes_Gauss_pla,1),1);
Jpla(ind_Gauss_INpla) = solk1.Jphi(ind_Gauss_INpla);

solk1.Jphi = Jpla;

% % F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),solk1.Jphi);
% % J_interp = zeros(meshData_loc.nn,1);
% % J_interp = F(meshData_loc.n(:,1),meshData_loc.n(:,2));
% % 
% % Jpla.faces=meshData_loc.t(meshData_loc.ind_t_InFW,1:3);
% % Jpla.vertices=meshData_loc.n;
% % Jpla.facevertexcdata=J_interp;
% % 
% % figure
% % edges=pdemesh(meshData.n',meshData.e_interface',[]);
% % set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
% % axis([Rplot Zplot]),  hold on,
% % hh=patch(Jpla,'facecolor','interp','edgecolor','none');
% % axis equal;  hold on; colorbar vert; % title('j_\phi');
% % axis([Rplot Zplot])
% % 
% % figure
% % plot3(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),solk1.Jphi,'o')

%%
fprintf(' \n');
fprintf('--> Done !!! \n');
time_end = toc(time_start_code);

fprintf('    TOTAL Eelapsed time is %5.1f seconds \n',  time_end)



%% Save output

if any(strcmp('GG_BEM',fieldnames(meshData_loc)))
    meshData_loc = rmfield(meshData_loc,'GG_BEM');
end

OUTPUT_FRIDA.plaparameter       = plaparameter;
% % OUTPUT_FRIDA.meshData           = meshData;
% % OUTPUT_FRIDA.meshData_loc       = meshData_loc;
% % OUTPUT_FRIDA.Conductors         = Conductors;
OUTPUT_FRIDA.solk1              = solk1;
OUTPUT_FRIDA.scarto_tot         = scarto_tot;
OUTPUT_FRIDA.residuo_tot        = residuo_tot;
OUTPUT_FRIDA.residuo_norm_tot   = residuo_norm_tot;
OUTPUT_FRIDA.time_per_iteration = time_per_iteration;
OUTPUT_FRIDA.time_solver        = time_solver;
OUTPUT_FRIDA.ii_freeB           = ii_freeB;
OUTPUT_FRIDA.ind_figure         = ind_figure;
% % OUTPUT_FRIDA.POSTPROC           = POSTPROC;


if SETTINGS.SAVE_OUTPUT == true
    
    save(['OUTPUT_FRIDA' SETTINGS.filename_out],'OUTPUT_FRIDA')
    
end

% % restoredefaultpath



