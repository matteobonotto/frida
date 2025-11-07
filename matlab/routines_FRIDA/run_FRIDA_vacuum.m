function OUT = run_FRIDA_vacuum(...
    meshData_ext, ...
    IE_evol, ...
    Conductors, ...
    SETTINGS)

n_time = length(IE_evol.time_sim);
% % plaparameter = IE_evol.plaparameter;

%% Output quantities


if ~any(strcmp('J_c_t0',fieldnames(SETTINGS)))
    SETTINGS.J_c_t0 = zeros(meshData_ext.n_pas,1);
end
if SETTINGS.IS_EVOL
    if ~any(strcmp('v_gap_t0',fieldnames(SETTINGS)))
        SETTINGS.v_gap_t0 = zeros(meshData_ext.n_cut_pas,1);
    end
end

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

if any(strcmp('flux_loops',fieldnames(SETTINGS)))
    Meas_fluxloops_TD = zeros(size(SETTINGS.flux_loops,1),n_time);
end

if any(strcmp('pickup',fieldnames(SETTINGS)))
    Meas_pickup_TD = zeros(size(SETTINGS.pickup,1),n_time);
end

if SETTINGS.IS_EVOL
    J_c_TD = zeros(meshData_ext.n_pas,n_time);
    v_gap_TD = zeros(meshData_ext.n_cut_pas,n_time);
end

% % Psi_TD = zeros(meshData_pla.nn,n_time);
% % I_pla_TD = zeros(meshData_pla.nn,n_time);
% % J_pla_TD = zeros(meshData_pla.n_Gauss*meshData_pla.nt,n_time);



%% First solution (t = t0)


%%% Define input parameters at t0
Conductors_t0 = Conductors;
Conductors_t0.Currents = Conductors_t0.Currents(:,1);

plaparameter_t0 = [];

% % plaparameter_t0.Centroid = plaparameter.Centroid;
% % plaparameter_t0.R_0      = plaparameter.R_0(1);
% % plaparameter_t0.beta_0   = plaparameter.beta_0(1);
% % plaparameter_t0.alpha_M  = plaparameter.alpha_M(1);
% % plaparameter_t0.alpha_N  = plaparameter.alpha_N(1);
% % plaparameter_t0.Ipla     = plaparameter.Ipla(1);


if ~any(strcmp('RR',fieldnames(SETTINGS)))
    Rplot=[min(meshData_ext.n(:,1))...
        max(meshData_ext.n(:,1))]+[-0.3 0.3];
    Zplot=[min(meshData_ext.n(:,2))...
        max(meshData_ext.n(:,2))]+[-0.3 0.3];

    rgrid=linspace(Rplot(1),Rplot(2),100);
    zgrid=linspace(Zplot(1),Zplot(2),100);

    [RR,ZZ] = meshgrid(rgrid,zgrid);

    SETTINGS.RR = RR;
    SETTINGS.ZZ = ZZ;
end
    
    
if SETTINGS.RUN == 0 
    %% VACUUM STATIC
    if any(strcmp('flux_loops',fieldnames(SETTINGS)))
        Meas_fluxloops_TD_t0 = meshData_ext.G_flux_VI_act*Conductors_t0.Currents;
        OUT.Meas_fluxloops = Meas_fluxloops_TD_t0;
    end
    
    if any(strcmp('pickup',fieldnames(SETTINGS)))
        Meas_pickup_TD_t0 = meshData_ext.G_B_pickup_VI_act*Conductors_t0.Currents;
        OUT.Meas_pickup = Meas_pickup_TD_t0;

    end
    
    tri_act = meshData_ext.t;
    nodes_act = meshData_ext.n;
    keyreg_act = meshData_ext.type;
    ind_act = meshData_ext.ind_act;
    n_act = length(ind_act);
    
    degree_G_source = 8;
    
    G_flux_VI_act =fun_G_flux_active_VI_stable_fast(n_act, ...
        tri_act, ...
        nodes_act, ...
        ind_act, ...
        keyreg_act, ...
        degree_G_source, ...
        RR(:), ...
        ZZ(:));
    
    OUT.Psi_RR = reshape(...
        G_flux_VI_act*Conductors_t0.Currents,...
        size(SETTINGS.RR,1),...
        size(SETTINGS.ZZ,1));
    
    [G_Br_VI_act,G_Bz_VI_act] =fun_G_BrBz_active_VI_stable_fast(n_act, ...
        tri_act, ...
        nodes_act, ...
        ind_act, ...
        keyreg_act, ...
        degree_G_source, ...
        RR(:), ...
        ZZ(:));
    
    OUT.Br_RR = reshape(...
        G_Br_VI_act*Conductors_t0.Currents,...
        size(SETTINGS.RR,1),...
        size(SETTINGS.ZZ,1));
    
    OUT.Br_ZZ = reshape(...
        G_Bz_VI_act*Conductors_t0.Currents,...
        size(SETTINGS.RR,1),...
        size(SETTINGS.ZZ,1));
    
    OUT.RR = SETTINGS.RR;
    OUT.ZZ = SETTINGS.ZZ;

else
    %% VACUUM EVOL
    uu = Conductors.Currents;
    time_sim = Conductors.time_sim;
    n_cut_pas = meshData_ext.n_cut_pas;
    
    x_t0 = -meshData_ext.A\meshData_ext.B_act*uu(:,1);
    tic
    xx_VI = fun_DAE_Crank_Nicolson_State_Space(...
        time_sim,...
        x_t0,...
        uu,...
        meshData_ext.E,...
        meshData_ext.A,...
        meshData_ext.B_act);
    toc

    a_VI = xx_VI(1:end-n_cut_pas,:);
    v_gap = xx_VI(end-n_cut_pas+1:end,:);

    jc = meshData_ext.C(:,1:end-n_cut_pas)*a_VI + meshData_ext.D_act*uu;
    
    % %     OUT.a_VI = a_VI;
    OUT.v_gap_TD = v_gap;
    OUT.J_c_TD = jc;
    
    %%% Compute measures (if any)
    if any(strcmp('flux_loops',fieldnames(SETTINGS)))
        Meas_fluxloops_TD = meshData_ext.G_flux_VI_act*Conductors_t0.Currents + ...
            meshData_ext.G_flux_VI_pas*jc;
        OUT.Meas_fluxloops_TD = Meas_fluxloops_TD;
    end
    
    if any(strcmp('pickup',fieldnames(SETTINGS)))
        Meas_pickup_TD = meshData_ext.G_B_pickup_VI_act*Conductors_t0.Currents + ...
            meshData_ext.G_B_pickup_VI_pas*jc;
        OUT.Meas_pickup_TD = Meas_pickup_TD;

    end
    
    
    %%% Compute Psi, Br,Bz on a structured grid
    tri_act = meshData_ext.t;
    nodes_act = meshData_ext.n;
    keyreg_act = meshData_ext.type;
    ind_act = meshData_ext.ind_act;
    n_act = length(ind_act);
    
    meshData_pas = fun_get_meshData_pas(meshData_ext);
    tri = meshData_pas.t;
    nodes = meshData_pas.n;
    
    degree_G_source = 4;
    nodes_matrix = [nodes(tri(:,1),:) ...
        nodes(tri(:,2),:) ...
        nodes(tri(:,3),:)];
    P1 = nodes_matrix(:,[1 2]);
    P2 = nodes_matrix(:,[3 4]);
    P3 = nodes_matrix(:,[5 6]);
    [~,P_G_source,n_G_source] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source);
    
    %%% Psi
    G_flux_VI_act =fun_G_flux_active_VI_stable_fast(n_act, ...
        tri_act, ...
        nodes_act, ...
        ind_act, ...
        keyreg_act, ...
        degree_G_source, ...
        RR(:), ...
        ZZ(:));
    
    G_flux_VI_pas = fun_G_flux_passive_VI_stable_fast(tri, ...
        nodes, ...
        N_order, ...
        degree_G_source, ...
        n_G_source, ...
        P_G_source,...
        RR(:), ...
        ZZ(:));
    
    OUT.Psi_RR = zeros(size(SETTINGS.RR,1),size(SETTINGS.ZZ,2),n_time);
    for i = 1:n_time
        OUT.Psi_RR(:,:,i) = reshape(...
        G_flux_VI_act*Conductors.Currents(:,i) + G_flux_VI_pas*OUT.J_c_TD(:,i),...
        size(SETTINGS.RR,1),...
        size(SETTINGS.ZZ,1));
    end
    
    %%% Br, Bz
    [G_Br_VI_pas,G_Bz_VI_pas] = fun_G_BrBz_passive_VI_stable_fast(tri, ...
        nodes, ...
        N_order, ...
        degree_G_source, ...
        n_G_source, ...
        P_G_source,...
        RR(:), ...
        ZZ(:));
    
    [G_Br_VI_act,G_Bz_VI_act] =fun_G_BrBz_active_VI_stable_fast(n_act, ...
        tri_act, ...
        nodes_act, ...
        ind_act, ...
        keyreg_act, ...
        degree_G_source, ...
        RR(:), ...
        ZZ(:));
    
    OUT.Br_RR = zeros(size(SETTINGS.RR,1),size(SETTINGS.ZZ,2),n_time);
    OUT.Br_ZZ = zeros(size(SETTINGS.RR,1),size(SETTINGS.ZZ,2),n_time);
    for i = 1:n_time
        OUT.Br_RR(:,:,i) = reshape(...
            G_Br_VI_act*Conductors.Currents(:,i) + G_Br_VI_pas*OUT.J_c_TD(:,i),...
            size(SETTINGS.RR,1),...
            size(SETTINGS.ZZ,1));
        OUT.Br_ZZ(:,:,i) = reshape(...
            G_Bz_VI_act*Conductors.Currents(:,i) + G_Bz_VI_pas*OUT.J_c_TD(:,i),...
            size(SETTINGS.RR,1),...
            size(SETTINGS.ZZ,1));
    end

    OUT.RR = SETTINGS.RR;
    OUT.ZZ = SETTINGS.ZZ;
    
end


























