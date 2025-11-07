function [OUT] = main_run_FRIDA_evol_VI(SETTINGS)

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

fprintf('   FRee-boundary Integro-Differential Axisymmetric (solver) - Time Domain\n\n')
fprintf('   M. Bonotto, D. Abate\n')
fprintf('   Consorzio RFX \n\n')
fprintf('   4/2024\n')
fprintf('   ver. 3.3\n')

fprintf('\n\n\n')
fprintf('\n\n\n')

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')



%% Some preliminary stuff
% VIP quantity
mu0=4*pi*1.e-7;


%%% set Latex as default figure font/interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');



%% Loading data and pre-processing

%%% Load evol data
% % load(['INPUT_FRIDA_evol_VI_' SETTINGS.filename_evol],'IE_evol')
% % Conductors   = IE_evol.Conductors;
% % plaparameter = IE_evol.plaparameter;

% load equil data
% % load(['./data_run_FRIDA/INPUT_FRIDA_equil_' SETTINGS.filename_equil], ...
load(['./data_in_FRIDA/INPUT_FRIDA_equil_', SETTINGS.filename, '.mat'], ...
    'IE_evol')
plaparameter = IE_evol.plaparameter;
Conductors = IE_evol.Conductors;
% % IE_evol.plaparameter = plaparameter;
% % IE_evol.Conductors = Conductors;
% % IE_evol.time_sim = Conductors.time_sim;

% load geometry data
% % load(['./data_run_FRIDA/INPUT_FRIDA_geo_' SETTINGS.filename_geo], ...
load(['./data_in_FRIDA/INPUT_FRIDA_geo_', SETTINGS.filename, '.mat'], ...
    'meshData_pla', ...
    'meshData', ...
    'sensors')
meshData_ext = meshData;
if isfield(sensors, 'pickup')
    SETTINGS.pickup = sensors.pickup;
end
if isfield(sensors, 'pickup')
    SETTINGS.flux_loops = sensors.flux_loops;
end



%% RUN TYPE
if SETTINGS.RUN == 2 || SETTINGS.RUN == 3 % evol vacuum or evol plasma
    SETTINGS.IS_EVOL = 1;
    n_time = length(IE_evol.time_sim);
else
    SETTINGS.IS_EVOL = 0;
    n_time = 1;
end


%%
%%% Initialization of SETTINGS and INPUT
run_initialize_SETTINGS


if SETTINGS.IS_EVOL
    time_sim     = plaparameter.time_sim;
end

%%% Pre-processing
if SETTINGS.PREPROC == true
    fprintf('\n\n\n')
    fprintf('==============================================\n')
    fprintf('*** RUNNING PRE-PROCESSING *** \n')
    
    time_start_preproc = tic;
    
    % load geometry
    %     load(['INPUT_FRIDA_evol_VI_geo_' SETTINGS.filename_geo], ...
    %         'meshData_pla', ...
    %         'meshData_ext')
    
    % run preprocessing
    run_FRIDA_evol_preprocessing
    
    fprintf(' \n');
    fprintf('*** --> PRE-PROCESSING DONE!!! *** \n');
    
    time_end_preproc = toc(time_start_preproc);
    
    fprintf('    TOTAL Eelapsed time is %5.1f seconds \n',  time_end_preproc)
    fprintf('==============================================\n')
    fprintf('\n\n')
    
else
    fprintf('\n\n\n')
    fprintf('==============================================\n')
    fprintf('*** SKIPPING PRE-PROCESSING *** \n')
    
    if SETTINGS.IS_EVOL
        load(['./data_in_FRIDA/INPUT_FRIDA_geo_preproc_', SETTINGS.filename, '.mat'], ...
            'meshData_pla', ...
            'meshData_ext',...
            'meshData_pas')
    else
        load(['./data_in_FRIDA/INPUT_FRIDA_geo_preproc_', SETTINGS.filename, '.mat'], ...
            'meshData_pla', ...
            'meshData_ext')
    end
    
    fprintf('\n')
    fprintf('==============================================\n')
    fprintf('\n\n')
    
end


%% Matrix for Crank-Nicolson

if SETTINGS.IS_EVOL
    fprintf('\n')
    disp('Computing matrices for Crank�Nicolson method ...')
    
    tic
    compute_Crank_Nicolson_matrices
    toc
end


%% Run simulation 

if SETTINGS.RUN == 0 || SETTINGS.RUN == 2 % Vacuum static || Evol
    OUT = run_FRIDA_vacuum(...
        meshData_ext, ...
        IE_evol, ...
        Conductors, ...
        SETTINGS);
    
elseif SETTINGS.RUN == 1 || SETTINGS.RUN == 3 % Plasma static || Evol
    OUT = run_FRIDA_plasma(...
        meshData_pla, ...
        meshData_ext, ...
        IE_evol, ...
        Conductors, ...
        SETTINGS);
    
end

% % 
% % %% Output quantities
% % 
% % 
% % if ~any(strcmp('J_c_t0',fieldnames(SETTINGS)))
% %     SETTINGS.J_c_t0 = zeros(meshData_ext.n_pas,1);
% % end
% % if SETTINGS.IS_EVOL
% %     if ~any(strcmp('v_gap_t0',fieldnames(SETTINGS)))
% %         SETTINGS.v_gap_t0 = zeros(meshData_ext.n_cut_pas,1);
% %     end
% % end
% % 
% % Separatrix_r_TD = zeros(200,n_time);
% % Separatrix_z_TD = zeros(200,n_time);
% % residuo_TD = zeros(1,n_time);
% % Centroid_TD = zeros(2,n_time);
% % Psi_B_TD = zeros(1,n_time);
% % Psi_Ax_TD = zeros(1,n_time);
% % Ax_r_TD = zeros(1,n_time);
% % Ax_z_TD = zeros(1,n_time);
% % Psi_Centroid_TD = zeros(1,n_time);
% % lambda_TD = zeros(1,n_time);
% % ii_freeB_TD = zeros(1,n_time);
% % 
% % if any(strcmp('flux_loops',fieldnames(SETTINGS)))
% %     Meas_fluxloops_TD = zeros(size(SETTINGS.flux_loops,1),n_time);
% % end
% % 
% % if any(strcmp('pickup',fieldnames(SETTINGS)))
% %     Meas_pickup_TD = zeros(size(SETTINGS.pickup,1),n_time);
% % end
% % 
% % if SETTINGS.IS_EVOL
% %     J_c_TD = zeros(meshData_ext.n_pas,n_time);
% %     v_gap_TD = zeros(meshData_ext.n_cut_pas,n_time);
% % end
% % 
% % Psi_TD = zeros(meshData_pla.nn,n_time);
% % I_pla_TD = zeros(meshData_pla.nn,n_time);
% % J_pla_TD = zeros(meshData_pla.n_Gauss*meshData_pla.nt,n_time);
% % 
% % 
% % 
% % %% First solution (t = t0)
% % 
% % 
% % %%% Define input parameters at t0
% % Conductors_t0 = Conductors;
% % Conductors_t0.Currents = Conductors_t0.Currents(:,1);
% % 
% % plaparameter_t0.Centroid = plaparameter.Centroid;
% % plaparameter_t0.R_0      = plaparameter.R_0(1);
% % plaparameter_t0.beta_0   = plaparameter.beta_0(1);
% % plaparameter_t0.alpha_M  = plaparameter.alpha_M(1);
% % plaparameter_t0.alpha_N  = plaparameter.alpha_N(1);
% % plaparameter_t0.Ipla     = plaparameter.Ipla(1);
% % 
% % 
% % if SETTINGS.ii_START == 1
% %     
% %     fprintf('\n\n')
% %     fprintf('==================================================================\n')
% %     fprintf('==================================================================\n')
% %     
% %     fprintf('                *** FIRST SOLUTION (t = t0) ***\n')
% %     
% %     fprintf('==================================================================\n')
% %     fprintf('==================================================================\n')
% %     fprintf('\n\n')
% %     
% %     %%% If passive current at t0 is not defined (SETTINGS.J_c_t0), set it to zero
% %     if SETTINGS.IS_EVOL
% %         J_c_TD(:,1) = SETTINGS.J_c_t0;
% %     end
% %     
% %     
% %     %%% Compute vacuum field by only BCs of active and passive current
% %     Psi_MS_bc_act = meshData_ext.G_BCs_VI_act*Conductors_t0.Currents;
% %     Psi_MS_bc_pas = meshData_ext.G_BCs_VI_pas*SETTINGS.J_c_t0;
% %     
% %     SETTINGS.Psi_MS_bc = (Psi_MS_bc_act + Psi_MS_bc_pas)/2/pi;
% %     SETTINGS.VAC_METHOD = 3;
% %     
% %     
% %     %%% Compute eqiulibrium at t0
% %     %%% STATIC VACUUM
% %     if SETTINGS.RUN == 0
% %         [OUTPUT_FRIDA_t0] = run_FRIDA_vacuum_t0(...
% %             meshData_pla,...
% %             meshData_ext,...
% %             Conductors_t0,...
% %             SETTINGS);
% %         
% %         Psi_TD(:,1) = 2*pi*solk_t0.Psi; % TOTAL poloidal flux
% %         I_pla_TD(:,1) = solk_t0.Iphi;
% %         J_pla_TD(:,1) = solk_t0.Jphi;
% %         Separatrix_r_TD(:,1) = solk_t0.Separatrix(:,1);
% %         Separatrix_z_TD(:,1) = solk_t0.Separatrix(:,2);
% %         residuo_TD(1) = OUTPUT_FRIDA_t0.solk1.residuo_norm;
% %         Centroid_TD(:,1) = [OUTPUT_FRIDA_t0.solk1.Centroid_RR OUTPUT_FRIDA_t0.solk1.Centroid_ZZ].';
% %         Psi_B_TD(:,1) = 2*pi*solk_t0.Psi_B; % TOTAL poloidal flux
% %         Psi_Ax_TD(:,1) = 2*pi*solk_t0.Psi_a; % TOTAL poloidal flux
% %         Psi_Centroid_TD(:,1) = 2*pi*solk_t0.Centroid_Psi; % TOTAL poloidal flux
% %         lambda_TD(:,1) = solk_t0.lambda;
% %         Ax_r_TD(:,1) = OUTPUT_FRIDA_t0.solk1.Axis_RR;
% %         Ax_z_TD(:,1) = OUTPUT_FRIDA_t0.solk1.Axis_ZZ;
% %         ii_freeB_TD(1) = OUTPUT_FRIDA_t0.ii_freeB;
% %         
% %         if any(strcmp('flux_loops',fieldnames(SETTINGS)))
% %             Meas_fluxloops_TD_t0 = meshData_ext.G_flux_VI_act*Conductors_t0.Currents + ...
% %                 meshData_pla.G_flux_loops_pla*OUTPUT_FRIDA_t0.solk1.Iphi + ...
% %                 meshData_ext.G_flux_VI_pas*SETTINGS.J_c_t0;
% %             
% %             Meas_fluxloops_TD(:,1) = Meas_fluxloops_TD_t0;
% %         end
% %         
% %         if any(strcmp('pickup',fieldnames(SETTINGS)))
% %             Meas_pickup_TD_t0 = meshData_ext.G_B_pickup_VI_act*Conductors_t0.Currents + ...
% %                 meshData_pla.G_pickup_pla*OUTPUT_FRIDA_t0.solk1.Iphi + ...
% %                 meshData_ext.G_B_pickup_VI_pas*SETTINGS.J_c_t0;
% %             Meas_pickup_TD(:,1) = Meas_pickup_TD_t0;
% %         end
% %         
% %     %%% STSTIC WITH PLASMA
% %     elseif SETTINGS.RUN == 1
% %         [OUTPUT_FRIDA_t0] = run_FRIDA_evol_t0(...
% %             meshData_pla,...
% %             meshData_ext,...
% %             Conductors_t0,...
% %             plaparameter_t0,...
% %             SETTINGS);
% %         
% %         % %     openfig kjsbdfskjbfds
% %         % %     plot(OUTPUT_FRIDA_t0.solk1.Separatrix(:,1),OUTPUT_FRIDA_t0.solk1.Separatrix(:,2),'-*k');
% %         
% %         %%% Post-processing
% %         solk_t0 = OUTPUT_FRIDA_t0.solk1;
% %         solk_t0.xx_th = solk_t0.xx_k;
% %         
% %         solk_t0.v_gap = SETTINGS.v_gap_t0;
% %         solk_t0.J_c_th = SETTINGS.J_c_t0;
% %         
% %         Psi_TD(:,1) = 2*pi*solk_t0.Psi; % TOTAL poloidal flux
% %         I_pla_TD(:,1) = solk_t0.Iphi;
% %         J_pla_TD(:,1) = solk_t0.Jphi;
% %         Separatrix_r_TD(:,1) = solk_t0.Separatrix(:,1);
% %         Separatrix_z_TD(:,1) = solk_t0.Separatrix(:,2);
% %         residuo_TD(1) = OUTPUT_FRIDA_t0.solk1.residuo_norm;
% %         Centroid_TD(:,1) = [OUTPUT_FRIDA_t0.solk1.Centroid_RR OUTPUT_FRIDA_t0.solk1.Centroid_ZZ].';
% %         Psi_B_TD(:,1) = 2*pi*solk_t0.Psi_B; % TOTAL poloidal flux
% %         Psi_Ax_TD(:,1) = 2*pi*solk_t0.Psi_a; % TOTAL poloidal flux
% %         Psi_Centroid_TD(:,1) = 2*pi*solk_t0.Centroid_Psi; % TOTAL poloidal flux
% %         lambda_TD(:,1) = solk_t0.lambda;
% %         Ax_r_TD(:,1) = OUTPUT_FRIDA_t0.solk1.Axis_RR;
% %         Ax_z_TD(:,1) = OUTPUT_FRIDA_t0.solk1.Axis_ZZ;
% %         ii_freeB_TD(1) = OUTPUT_FRIDA_t0.ii_freeB;
% %         
% %         if any(strcmp('flux_loops',fieldnames(SETTINGS)))
% %             Meas_fluxloops_TD_t0 = meshData_ext.G_flux_VI_act*Conductors_t0.Currents + ...
% %                 meshData_pla.G_flux_loops_pla*OUTPUT_FRIDA_t0.solk1.Iphi + ...
% %                 meshData_ext.G_flux_VI_pas*SETTINGS.J_c_t0;
% %             
% %             Meas_fluxloops_TD(:,1) = Meas_fluxloops_TD_t0;
% %         end
% %         
% %         if any(strcmp('pickup',fieldnames(SETTINGS)))
% %             Meas_pickup_TD_t0 = meshData_ext.G_B_pickup_VI_act*Conductors_t0.Currents + ...
% %                 meshData_pla.G_pickup_pla*OUTPUT_FRIDA_t0.solk1.Iphi + ...
% %                 meshData_ext.G_B_pickup_VI_pas*SETTINGS.J_c_t0;
% %             Meas_pickup_TD(:,1) = Meas_pickup_TD_t0;
% %         end
% %     end
% %     
% % end
% % 
% % %% Time-stepping scheme
% % 
% % if SETTINGS.IS_EVOL
% %     
% %     fprintf('\n\n')
% %     fprintf('==================================================================\n')
% %     fprintf('==================================================================\n')
% %     
% %     % % fprintf('\n')
% %     fprintf('             *** STARTING TIME DOMAIN SIMULATION ***\n')
% %     % % fprintf('\n\n')
% %     
% %     fprintf('==================================================================\n')
% %     fprintf('==================================================================\n')
% %     fprintf('\n\n')
% %     
% %     SETTINGS.RUN_FIXBOUNDARY = false;
% %     
% %     v_gap_TD(:,1) = SETTINGS.v_gap_t0;
% %     
% %     
% %     if SETTINGS.ii_START > 1
% %         filename_load = sprintf('./temp_out/input_ii_%i.mat',SETTINGS.ii_START);
% %         load(filename_load, 'solk_old')
% %     else
% %         solk_old = solk_t0;
% %     end
% %     
% %     
% %     lambda_Ipla_th = solk_old.Iphi(meshData_pla.ind_D);
% %     a_pas_t0 = L_VI*solk_old.J_c_th + ...
% %         M_act_pas*Conductors_t0.Currents + ...
% %         M_CS_pas*T_Ipla_eq_D*lambda_Ipla_th;
% %     
% %     xx_DAE_t0 = [a_pas_t0; solk_old.v_gap];
% %     solk_old.xx_DAE_old = xx_DAE_t0;
% %     
% %     J_c_th = solk_old.J_c_th;
% %     
% %     if isfield(solk_old,'xx_k'); solk_old.xx_th = solk_old.xx_k; solk_old = rmfield(solk_old,'xx_k'); end
% %     
% %     SETTINGS.VAC_METHOD = 3;
% %     SETTINGS.FIRST_SOLUTION = true;
% %     SETTINGS.FIRST_SOLUTION = true;
% %     SETTINGS.NFREEB_FIXB = 0;
% %     SETTINGS.TOLL = 1e-10;
% %     SETTINGS.SCARTO =  1e-10;
% %     SETTINGS.SOLVER_RELAXED = true;
% %     % % SETTINGS.FIGURES = true;
% %     SETTINGS.ii_freeB_max = 25;
% %     SETTINGS.SAVE_OUTPUT = false;
% %     SETTINGS.SOLVER = 'NR';
% %     SETTINGS.ii_CN_tot = n_time;
% %     
% %     
% %     %%%
% %     time_CN_start = tic;
% %     hh_CN = waitbar(0,sprintf('Running Crank-Nicolson: step %i of %i',SETTINGS.ii_START+1,n_time));
% %     for ii_CN = SETTINGS.ii_START+1:n_time
% %         
% %         waitbar(ii_CN/n_time,hh_CN,sprintf('Running Crank-Nicolson: step %i of %i',ii_CN,n_time))
% %         SETTINGS.ii_CN = ii_CN;
% %         
% %         ind_time = ii_CN;
% %         
% %         plaparameter_th.Centroid = [2.01 -.04];
% %         plaparameter_th.Ipla     = IE_evol.plaparameter.Ipla(ind_time);
% %         plaparameter_th.R_0      = IE_evol.plaparameter.R_0(ind_time);
% %         plaparameter_th.beta_0   = IE_evol.plaparameter.beta_0(ind_time);
% %         plaparameter_th.alpha_M  = IE_evol.plaparameter.alpha_M(ind_time);
% %         plaparameter_th.alpha_N  = IE_evol.plaparameter.alpha_N(ind_time);
% %         
% %         Conductors_th.Nconductors  = IE_evol.Conductors.Nconductors(1);
% %         Conductors_th.Nturns       = IE_evol.Conductors.Nturns(:,1);
% %         Conductors_th.Currents     = IE_evol.Conductors.Currents(:,ind_time);%.*IE_evol.Conductors.Nturns(:,2);
% %         Conductors_th.Currents_old = IE_evol.Conductors.Currents(:,ind_time-1);%.*IE_evol.Conductors.Nturns(:,2);
% %         
% %         %%% Compute vacuum field by only BCs of active and passive current
% %         Psi_MS_bc_act = meshData_ext.G_BCs_VI_act*Conductors_th.Currents;
% %         Psi_MS_bc_pas = meshData_ext.G_BCs_VI_pas*solk_old.J_c_th;
% %         
% %         SETTINGS.Psi_MS_bc = (Psi_MS_bc_act + Psi_MS_bc_pas)/2/pi;
% %         SETTINGS.Psi_MS_bc_act = Psi_MS_bc_act/2/pi;
% %         
% %         SETTINGS.xx_0 = solk_old.xx_th;
% %         
% %         close all
% %         
% %         
% %         try
% % 
% %             solk_old_save = solk_old;
% %             if isfield(solk_old,'LU')
% %                 solk_old = rmfield(solk_old,'LU');
% %             end
% %             filename_save = sprintf('./temp_out/input_ii_%i.mat',ii_CN-1);
% %             solk_old = solk_old_save;
% % 
% %             [OUT_th] = run_FRIDA_TD_CN(meshData_pla,... % CN = Crank Nicolson
% %                 meshData_ext,...
% %                 plaparameter_th,...
% %                 Conductors_th,...
% %                 solk_old, ...
% %                 SETTINGS);
% %             
% %             if  OUT_th.ii_freeB >= SETTINGS.ii_freeB_max && ...
% %                     OUT_th.residuo_norm_tot(end) >= 5e-5
% %                 error('better to stop here!')
% %             end
% %             
% %             %%%
% %             solk_old.xx_th = OUT_th.solk1.xx_k;
% %             solk_old.J_c_th = OUT_th.solk1.J_c;
% %             
% %             %%%
% %             solk_old.xx_th           = OUT_th.solk1.xx_k;
% %             solk_old.xx_DAE_old      = OUT_th.solk1.xx_DAE;
% %             solk_old.Iphi            = OUT_th.solk1.Iphi;
% %             solk_old.Grad_nodes      = OUT_th.solk1.Grad_nodes;
% %             solk_old.Psi_Gauss       = OUT_th.solk1.Psi_Gauss;
% %             solk_old.Grad_Gauss      = OUT_th.solk1.Grad_Gauss;
% %             solk_old.Psi_c_t         = OUT_th.solk1.Psi_c_t;
% %             solk_old.Grad_c_t        = OUT_th.solk1.Grad_c_t;
% %             solk_old.Psi_axis        = OUT_th.solk1.Psi_axis;
% %             solk_old.Axis_RR         = OUT_th.solk1.Axis_RR;
% %             solk_old.Axis_ZZ         = OUT_th.solk1.Axis_ZZ;
% %             solk_old.ind_n_axis      = OUT_th.solk1.ind_n_axis;
% %             solk_old.delta_a_weights = OUT_th.solk1.delta_a_weights;
% %             solk_old.Psi_B           = OUT_th.solk1.Psi_B;
% %             solk_old.XP_RR           = OUT_th.solk1.XP_RR;
% %             solk_old.XP_ZZ           = OUT_th.solk1.XP_ZZ;
% %             solk_old.ind_t_Xp        = OUT_th.solk1.ind_t_Xp;
% %             solk_old.ind_n_XP        = OUT_th.solk1.ind_n_XP;
% %             solk_old.delta_b_weights = OUT_th.solk1.delta_b_weights;
% %             solk_old.Separatrix      = OUT_th.solk1.Separatrix;
% %             solk_old.Psi_Bpla        = OUT_th.solk1.Psi_Bpla;
% %             
% %             %%% For Quasi Newton
% %             if SETTINGS.QuasiNewton == true
% %                 solk_old.LU = OUT_th.solk1.LU;
% %                 solk_old.RR_psi_k = OUT_th.solk1.RR_psi_k;
% %             end
% %             
% %             %%% Measures
% %             if any(strcmp('flux_loops',fieldnames(SETTINGS)))
% %                 Meas_fluxloops_TD_th = meshData_ext.G_flux_VI_act*Conductors_th.Currents + ...
% %                     meshData_pla.G_flux_loops_pla*OUT_th.solk1.Iphi + ...
% %                     meshData_ext.G_flux_VI_pas*solk_old.J_c_th;
% %                 
% %                 Meas_fluxloops_TD(:,ii_CN) = Meas_fluxloops_TD_th;
% %             end
% %             
% %             if any(strcmp('pickup',fieldnames(SETTINGS)))
% %                 Meas_pickup_TD_th = meshData_ext.G_B_pickup_VI_act*Conductors_th.Currents + ...
% %                     meshData_pla.G_pickup_pla*OUT_th.solk1.Iphi + ...
% %                     meshData_ext.G_B_pickup_VI_pas*solk_old.J_c_th;
% %                 
% %                 Meas_pickup_TD(:,ii_CN) = Meas_pickup_TD_th;
% %             end
% %             
% %             
% %             %%% output quantities
% %             Separatrix_r_TD(:,ii_CN) = OUT_th.solk1.Separatrix(:,1);
% %             Separatrix_z_TD(:,ii_CN) = OUT_th.solk1.Separatrix(:,2);
% %             residuo_TD(ii_CN)        = OUT_th.solk1.residuo_norm;
% %             Centroid_TD(:,ii_CN)     = [OUT_th.solk1.Centroid_RR OUT_th.solk1.Centroid_ZZ].';
% %             J_c_TD(:,ii_CN)          = solk_old.J_c_th;
% %             v_gap_TD(:,ii_CN)        = OUT_th.solk1.v_gap;
% %             Psi_TD(:,ii_CN)          = 2*pi*OUT_th.solk1.Psi; % TOTAL poloidal flux
% %             Psi_B_TD(:,ii_CN)        = 2*pi*OUT_th.solk1.Psi_B; % TOTAL poloidal flux
% %             Psi_Ax_TD(:,ii_CN)       = 2*pi*OUT_th.solk1.Psi_a; % TOTAL poloidal flux
% %             Psi_Centroid_TD(:,ii_CN) = 2*pi*OUT_th.solk1.Centroid_Psi; % TOTAL poloidal flux
% %             I_pla_TD(:,ii_CN)        = OUT_th.solk1.Iphi;
% %             J_pla_TD(:,ii_CN)        = OUT_th.solk1.Jphi;
% %             lambda_TD(:,ii_CN)       = OUT_th.solk1.lambda;
% %             Ax_r_TD(:,ii_CN)         = OUT_th.solk1.Axis_RR;
% %             Ax_z_TD(:,ii_CN)         = OUT_th.solk1.Axis_ZZ;
% %             ii_freeB_TD(:,ii_CN)     = OUT_th.ii_freeB;
% %             
% %             
% %         catch me
% %             warning('Something went wrong. %s\n',me.message)
% %             fprintf('%s',me.message)
% %             % %             save data_debug_2
% %             aa = 0;
% %             % %         warning(me.identifier,'%s',me.message)
% %             break
% %             
% %         end
% %         
% %     end
% %     
% %     close(hh_CN)
% %     time_CN = toc(time_CN_start);
% %     
% % end
% % 
% % %% Output quantity
% % 
% % if ~SETTINGS.IS_EVOL
% %     time_CN = [];
% %     time_sim = 0;
% % end
% % 
% % if SETTINGS.ii_START > 1
% %     
% %     for ii_CN = 1:SETTINGS.ii_START
% %         filename_load = sprintf('./temp_out/input_ii_%i.mat',ii_CN);
% %         load(filename_load, 'solk_old')
% %         
% %         
% %         Separatrix_r_TD(:,ii_CN) = solk_old.Separatrix(:,1);
% %         Separatrix_z_TD(:,ii_CN) = solk_old.Separatrix(:,2);
% %         residuo_TD(ii_CN)        = solk_old.residuo_norm;
% %         Centroid_TD(:,ii_CN)     = [solk_old.Centroid_RR solk_old.Centroid_ZZ].';
% %         J_c_TD(:,ii_CN)          = solk_old.J_c_th;
% %         v_gap_TD(:,ii_CN)        = solk_old.v_gap;
% %         Psi_TD(:,ii_CN)          = 2*pi*solk_old.Psi; % TOTAL poloidal flux
% %         Psi_B_TD(:,ii_CN)        = 2*pi*solk_old.Psi_B; % TOTAL poloidal flux
% %         Psi_Ax_TD(:,ii_CN)       = 2*pi*solk_old.Psi_a; % TOTAL poloidal flux
% %         Psi_Centroid_TD(:,ii_CN) = 2*pi*solk_old.Centroid_Psi; % TOTAL poloidal flux
% %         I_pla_TD(:,ii_CN)        = solk_old.Iphi;
% %         J_pla_TD(:,ii_CN)        = solk_old.Jphi;
% %         lambda_TD(:,ii_CN)       = solk_old.lambda;
% %         Ax_r_TD(:,ii_CN)         = solk_old.Axis_RR;
% %         Ax_z_TD(:,ii_CN)         = solk_old.Axis_ZZ;
% %         
% %         
% %         Conductors_th.Currents     = IE_evol.Conductors.Currents(:,ii);%.*IE_evol.Conductors.Nturns(:,2);
% %         
% %         if any(strcmp('flux_loops',fieldnames(SETTINGS)))
% %             Meas_fluxloops_TD_th = meshData_ext.G_flux_VI_act*Conductors_th.Currents + ...
% %                 meshData_pla.G_flux_loops_pla*solk_old.Iphi + ...
% %                 meshData_ext.G_flux_VI_pas*solk_old.J_c_th;
% %             
% %             Meas_fluxloops_TD(:,ii_CN) = Meas_fluxloops_TD_th;
% %         end
% %         
% %         if any(strcmp('pickup',fieldnames(SETTINGS)))
% %             Meas_pickup_TD_th = meshData_ext.G_B_pickup_VI_act*Conductors_th.Currents + ...
% %                 meshData_pla.G_pickup_pla*solk_old.Iphi + ...
% %                 meshData_ext.G_B_pickup_VI_pas*solk_old.J_c_th;
% %             
% %             Meas_pickup_TD(:,ii_CN) = Meas_pickup_TD_th;
% %         end
% %         
% %     end
% % end
% % 
% % 
% % 
% % %%%
% % n_time_new = numel(find(Ax_r_TD ~= 0));
% % ii_time_new = 1:n_time_new;
% % 
% % OUT.time_sim          = time_sim(ii_time_new);
% % OUT.Conductors        = Conductors;
% % OUT.IE_evol           = IE_evol;
% % 
% % OUT.Separatrix_r_TD   = Separatrix_r_TD(:,ii_time_new);
% % OUT.Separatrix_z_TD   = Separatrix_z_TD(:,ii_time_new);
% % OUT.residuo_TD        = residuo_TD(ii_time_new);
% % OUT.Centroid_TD       = Centroid_TD(:,ii_time_new);
% % OUT.J_c_TD            = J_c_TD(:,ii_time_new);
% % OUT.v_gap_TD          = v_gap_TD(:,ii_time_new);
% % OUT.Psi_TD            = Psi_TD(:,ii_time_new);
% % OUT.Psi_B_TD          = Psi_B_TD(ii_time_new);
% % OUT.Psi_Ax_TD         = Psi_Ax_TD(ii_time_new);
% % OUT.Ax_r_TD           = Ax_r_TD(ii_time_new);
% % OUT.Ax_z_TD           = Ax_z_TD(ii_time_new);
% % OUT.Psi_Centroid_TD   = Psi_Centroid_TD(ii_time_new);
% % OUT.I_pla_TD          = I_pla_TD(:,ii_time_new);
% % OUT.J_pla_TD          = J_pla_TD(:,ii_time_new);
% % OUT.lambda_TD         = lambda_TD(ii_time_new);
% % OUT.Meas_fluxloops_TD = Meas_fluxloops_TD(:,ii_time_new);
% % OUT.Meas_pickup_TD    = Meas_pickup_TD(:,ii_time_new);
% % 
% % OUT.time_CN           = time_CN;
% % OUT.ii_freeB_TD       = ii_freeB_TD;
% % 


end
































