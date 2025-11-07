

%

N_order = meshData_pla.shape_functions_N_order;

tri   = meshData_pla.t;
nodes = meshData_pla.n;
nn    = meshData_pla.nn;

%%
if ~any(strcmp('e_interface',fieldnames(meshData_ext)))
    
    fprintf('\n')
    disp('Computing boundary edge (for plotting purposes only) ...')
    
    ind_t_chosen = fun_vec_find(meshData_ext.type,meshData_ext.ind_act);
    t_chosen = meshData_ext.t(ind_t_chosen,:);
    % %     t_choosen=meshData_ext.t((sum(meshData_ext.type.' == meshData_ext.ind_act) ~= 0).',:);
    tic;
    if SETTINGS.RUN_MEX_ROUTINE
        [~,lati_bordo_act]=fun_lati_fast_mex(t_chosen,1);
    else
        [~,lati_bordo_act]=fun_lati_fast(t_chosen,1);
    end
    
    ind_t_chosen = fun_vec_find(meshData_ext.type,meshData_ext.ind_pas);
    t_chosen = meshData_ext.t(ind_t_chosen,:);
    % %     t_chosen=meshData_ext.t((sum(meshData_ext.type.' == meshData_ext.ind_pas) ~= 0).',:);
    if SETTINGS.RUN_MEX_ROUTINE
        [~,lati_bordo_pas]=fun_lati_fast_mex(t_chosen,1);
    else
        [~,lati_bordo_pas]=fun_lati_fast(t_chosen,1);
    end
    
    toc
    
    meshData_ext.e_interface = unique([lati_bordo_act; lati_bordo_pas],'rows');
    
end

if ~any(strcmp('e_interface',fieldnames(meshData_pla)))
    
    t_chosen = meshData_pla.t(meshData_pla.type >= 0,:);
    tic;
    [~,lati_bordo]=fun_lati_fast_mex(t_chosen,1);
    toc
    
    meshData_pla.e_interface = lati_bordo;
    
end


%% Domain and boundary nodes
meshData_pla.ind_B = meshData_pla.ind_n_bc;
meshData_pla.ind_D = setdiff(1:nn,meshData_pla.ind_B).';

meshData_pla.nn_B = length(meshData_pla.ind_B);
meshData_pla.nn_D = length(meshData_pla.ind_D);



%%

if SETTINGS.FIGURES == true
    figure
    hold on;    grid on
    
    mesh_vac = triplot(meshData_ext.t(:,1:3),meshData_ext.n(:,1),meshData_ext.n(:,2));
    set(mesh_vac,'color',.7*[1 1 1]);
    
    for ii = 1:Conductors.Nconductors
        mesh_vac=triplot(meshData_ext.t(meshData_ext.type == ii,1:3),meshData_ext.n(:,1),meshData_ext.n(:,2));
        set(mesh_vac,'color','r');
    end
    
    mesh_vac = triplot(meshData_pla.t(:,1:3),meshData_pla.n(:,1),meshData_pla.n(:,2));
    set(mesh_vac,'color',[30 75 93]*1e-2);
    
    edges=pdemesh(meshData_pla.n',meshData_pla.e_interface',[]);
    set(edges,'color','g','LineWidth',1);  axis equal;
    axis equal;  hold on;  box on; xlabel('r [m]'), ylabel('z [m]')
    
    edges=pdemesh(meshData_ext.n',meshData_ext.e_interface',[]);
    set(edges,'color','k')
    
    drawnow
    
end



%% Shape functions (for Integral and FEM)
tic

fprintf('\n')
disp('Computing shape functions ...')

if SETTINGS.RUN_MEX_ROUTINE == false
    [shape_f,shape_f_norm] = fun_shape_functions_stable(tri,nodes,N_order);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [shape_f,shape_f_norm] = fun_shape_functions_stable_mex(tri,nodes,N_order);
end

meshData_pla.shape_functions = shape_f;
meshData_pla.shape_functions_norm = shape_f_norm;

toc



%% Gauss nodes for plasma and stiffness matrix computation
fprintf('\n')
disp('Plasma domain: computing Gauss points for plasma (loc)...')
tic

% % n_Gauss = 3*N_order;
GAUSS_QUAD_DEGREE_PLA = SETTINGS.GAUSS_QUAD_DEGREE_PLA;

tri = meshData_pla.t;

nodes_Pla = [meshData_pla.n(tri(:,1),:) ...
    meshData_pla.n(tri(:,2),:) ...
    meshData_pla.n(tri(:,3),:)];

P1 = nodes_Pla(:,[1 2]);
P2 = nodes_Pla(:,[3 4]);
P3 = nodes_Pla(:,[5 6]);

if SETTINGS.RUN_MEX_ROUTINE == false
    [ww_Gauss_pla,nodes_Gauss_pla,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,GAUSS_QUAD_DEGREE_PLA);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [ww_Gauss_pla,nodes_Gauss_pla,n_Gauss] = fun_Gauss_points_triangle_stable_mex(P1,P2,P3,GAUSS_QUAD_DEGREE_PLA);
end

meshData_pla.n_Gauss = n_Gauss;
meshData_pla.nodes_Gauss_pla = nodes_Gauss_pla;
meshData_pla.ww_Gauss_pla = ww_Gauss_pla;



%% Incidence Matrix between nodes and triangles
fprintf('\n')
disp('Building incidence Matrix between nodes and triangles (loc)...')

tic

if SETTINGS.RUN_MEX_ROUTINE == false
    [MatInd_nodes_tri] = fun_MatInd_nodes_tri(N_order,tri,nn);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [MatInd_nodes_tri] = fun_MatInd_nodes_tri_fast_mex(N_order,tri,nn);
end

meshData_pla.MatInd_nodes_tri=MatInd_nodes_tri;

toc



%% Incidence Matrix between triangles and Gauss points
fprintf('\n')
disp('Building incidence Matrix between triangles and Gauss points (loc)...')
tic

index_cols = 1:meshData_pla.nt*meshData_pla.n_Gauss;
index_rows = reshape(repmat(1:meshData_pla.nt,meshData_pla.n_Gauss,1),meshData_pla.nt*meshData_pla.n_Gauss,1);

MatInd_tri_Gauss = reshape(index_cols,meshData_pla.n_Gauss,meshData_pla.nt)';

meshData_pla.MatInd_tri_Gauss=MatInd_tri_Gauss;

toc



%% Incidence Matrix between nodes
fprintf('\n')
disp('Building Incidence Matrix between nodes (loc)...')
tic

[MatInd_nodes] = fun_MatInd_nodes_v2(tri,nodes,MatInd_nodes_tri);
meshData_pla.MatInd_nodes=MatInd_nodes;

toc



%% Stiffness Matrix (FEM)
fprintf('\n')
disp('Computing stiffness matrix (loc)...')
tic

meshData_pla.GAUSS_QUAD_DEGREE_STIFFMAT = SETTINGS.GAUSS_QUAD_DEGREE_STIFFMAT;

% Already with MEX options!!!
[KK_nobc] = fun_StiffMat_FEM(meshData_pla,SETTINGS.RUN_MEX_ROUTINE);

meshData_pla.KK_nobc = KK_nobc;
toc



%% Boundary-Stiffness matrix (BEM)
fprintf('\n')
disp('Computing Green Matrix for BEM (pla -> CS) ...')

tic

clear source point
point.RR=meshData_pla.n(meshData_pla.ind_B,1);
point.ZZ=meshData_pla.n(meshData_pla.ind_B,2);
source.R=meshData_pla.n(meshData_pla.ind_D,1);
source.Z=meshData_pla.n(meshData_pla.ind_D,2);

if SETTINGS.RUN_MEX_ROUTINE == false
    
    source.current=ones(numel(source.R),1);
    source.type='ring';
    GG_BEM = fun_GreenMat_Loop_flux(source,point);
    GG_BEM.psi = GG_BEM.psi/2/pi;
    
elseif SETTINGS.RUN_MEX_ROUTINE == true
    
    R_points = point.RR(:);
    Z_points = point.ZZ(:);
    npt_point = length(R_points);
    
    OPT_PARALLEL = 0;
    n_threads = 24;
    
    GG_BEM.psi = zeros(npt_point,length(source.R));
    for ii=1:length(source.R)
        
        R_source = source.R(ii);
        Z_source = source.Z(ii);
        npt_source = 1;
        I_source = 1;
        
        psi = fun_Green_filament_flux_SP_f90(npt_source, ...
            R_source, ...
            Z_source, ...
            I_source, ...
            npt_point, ...
            R_points, ...
            Z_points, ...
            OPT_PARALLEL, ...
            n_threads);
        psi = psi/2/pi;
        
        GG_BEM.psi(:,ii) = psi;
        
    end
end

meshData_pla.GG_BEM = GG_BEM;
toc


%% Green matrix plasma -> probes
fprintf('\n')
disp('Computing Green matrix plasma -> probes ...')
OPT_PARALLEL = false;

tic
if any(strcmp('flux_loops',fieldnames(SETTINGS)))
    G_flux_loops_pla = fun_Green_pla_fil_flux_loops(meshData_pla,SETTINGS.flux_loops,SETTINGS,OPT_PARALLEL);
end

if any(strcmp('pickup',fieldnames(SETTINGS)))
    G_pickup_pla = fun_Green_pla_fil_pickup(meshData_pla,SETTINGS.pickup,SETTINGS,OPT_PARALLEL);
end
toc

meshData_pla.G_flux_loops_pla = G_flux_loops_pla;
meshData_pla.G_pickup_pla = G_pickup_pla;

%% Green matrix coils -> BCs nodes
% % fprintf('\n')
% % disp('Computing Green Matrix for BEM (coils -> CS) ...')
% %
% % n_CS = size(meshData_pla.C_surf,1);
% %
% % tic
% % Points_bc = meshData_pla.n(meshData_pla.ind_n_bc,:);
% % GG_tilde_coils_CS = zeros(n_CS,Conductors.Nconductors);
% %
% % for ii = 1:Conductors.Nconductors
% %
% %     ind_t_ii = find(meshData_ext.type == meshData_ext.ind_act(ii));
% %
% %     %
% %     ind_t_jj = ind_t_ii(:);
% %     tri_jj = meshData_ext.t(ind_t_jj,:);
% %
% %     nodes_jj = [meshData_ext.n(tri_jj(:,1),:) ...
% %         meshData_ext.n(tri_jj(:,2),:) ...
% %         meshData_ext.n(tri_jj(:,3),:)];
% %
% %     P1 = nodes_jj(:,[1 2]);
% %     P2 = nodes_jj(:,[3 4]);
% %     P3 = nodes_jj(:,[5 6]);
% %
% %     if SETTINGS.RUN_MEX_ROUTINE == false
% %         [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,SETTINGS.GAUSS_QUAD_DEGREE_PLA);
% %     elseif SETTINGS.RUN_MEX_ROUTINE == true
% %         [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable_mex(P1,P2,P3,SETTINGS.GAUSS_QUAD_DEGREE_PLA);
% %     end
% %
% %     %
% %     Ic = 1;
% %     Area_coil_ii = sum(ww);
% %     J_cond = Ic/Area_coil_ii;
% %     j_source = J_cond*ww;
% %
% %     clear source point
% %     point = Points_bc;
% %     source = nodes;
% %
% %     if SETTINGS.RUN_MEX_ROUTINE == false
% %
% %         Psi_MS_bc_ii = fun_Green_Flux_Loop( source, point, j_source);
% %         Psi_MS_bc_ii = Psi_MS_bc_ii/2/pi;
% %
% %     elseif SETTINGS.RUN_MEX_ROUTINE == true
% %
% %         R_points = point(:,1);
% %         Z_points = point(:,2);
% %         npt_point = length(R_points);
% %
% %         OPT_PARALLEL = 1;
% %         n_threads = 24;
% %
% %         R_source = source(:,1);
% %         Z_source = source(:,2);
% %         npt_source = length(R_source);
% %         I_source = j_source;
% %
% %         psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
% %         Psi_MS_bc_ii = psi/2/pi;
% %
% %     end
% %
% %     GG_tilde_coils_CS(:,ii) = Psi_MS_bc_ii;
% %
% % end
% % toc
% %
% % meshData_ext.GG_tilde_coils_CS = GG_tilde_coils_CS;


%% Matrix T_eq for plasma response
if SETTINGS.IS_EVOL == 1
    fprintf('\n')
    disp('Computing matrix T_eq for plasma response ...')
    
    C_surf = meshData_pla.C_surf;
    n_pla_eq = size(C_surf,1);
    
    centro = sum(C_surf)/size(C_surf,1);
    R0 = centro(1);
    Z0 = centro(2);
    
    [theta_C,rho_C] = cart2pol(C_surf(:,1)-R0,C_surf(:,2)-Z0);
    
    % define a target surface where the flux by I_p has to match the reference one
    r_target = 1.01*max(rho_C);
    C_target = [R0+r_target*cos(theta_C) r_target*sin(theta_C)];
    
    
    %%%
    tic
    if SETTINGS.RUN_MEX_ROUTINE == false
        
        clear source point
        point.RR=C_target(:,1);
        point.ZZ=C_target(:,2);
        source.R=C_surf(:,1);
        source.Z=C_surf(:,2);
        source.current=ones(numel(source.R),1);
        temp  = fun_GreenMat_Loop_flux( source, point );
        GG_eq_target = temp.psi;
        
        clear source point
        point.RR=C_target(:,1);
        point.ZZ=C_target(:,2);
        source.R=meshData_pla.n(:,1);
        source.Z=meshData_pla.n(:,2);
        source.current=ones(numel(source.R),1);
        temp = fun_GreenMat_Loop_flux(source,point);
        GG_pla_target = temp.psi;
        
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        
        %%% plasma nodes -> target surf
        point = C_target;
        source = meshData_pla.n;
        
        GG_pla_target = zeros(n_pla_eq,meshData_pla.nn);
        for ii = 1:meshData_pla.nn
            R_points = point(:,1);
            Z_points = point(:,2);
            npt_point = length(R_points);
            
            if n_pla_eq > 250
                OPT_PARALLEL = 1;
            else
                OPT_PARALLEL = 0;
            end
            n_threads = 24;
            
            R_source = source(ii,1);
            Z_source = source(ii,2);
            npt_source = 1;
            I_source = 1;
            
            psi_ii = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
            
            GG_pla_target(:,ii) = psi_ii;
        end
        
        
        %%% plasma CS -> target surf
        point = C_target;
        source = C_surf;
        
        GG_eq_target = zeros(n_pla_eq,n_pla_eq);
        for ii = 1:n_pla_eq
            R_points = point(:,1);
            Z_points = point(:,2);
            npt_point = length(R_points);
            
            if n_pla_eq > 250
                OPT_PARALLEL = 1;
            else
                OPT_PARALLEL = 0;
            end
            n_threads = 24;
            
            R_source = source(ii,1);
            Z_source = source(ii,2);
            npt_source = 1;
            I_source = 1;
            
            psi_ii = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
            
            GG_eq_target(:,ii) = psi_ii;
        end
        
    end
    toc
    
    T_Ipla_eq = GG_eq_target\GG_pla_target;
    % % T_Ipla_eq_D = T_Ipla_eq(:,meshData_pla.ind_D);
    
    meshData_pla.T_Ipla_eq = T_Ipla_eq;
end

%% VI matrices
if SETTINGS.IS_EVOL
    % % if SETTINGS.IS_EVOL
    
    fprintf('\n')
    disp('Computing VI matrices ...')
    
    [meshData_pas,...
        L_VI, ...
        R_VI, ...
        U_VI, ...
        V_VI, ...
        D_VI, ...
        M_VI, ...
        G_BCs_VI_act, ...
        G_BCs_VI_pas, ...
        G_CS_VI_act, ...
        G_CS_VI_pas, ...
        G_flux_VI_pas, ...
        G_B_pickup_VI_pas, ...
        G_flux_VI_act, ...
        G_B_pickup_VI_act, ...
        vec_A] = fun_compute_VI_matrices_FRIDA_evol(meshData_ext,meshData_pla,SETTINGS);
    
    n_pas = size(R_VI,1);
    meshData_ext.n_pas = n_pas;
    
    meshData_ext.L_VI              = L_VI;
    meshData_ext.R_VI              = R_VI;
    meshData_ext.U_VI              = U_VI;
    meshData_ext.V_VI              = V_VI;
    meshData_ext.D_VI              = D_VI;
    meshData_ext.M_VI              = M_VI;
    meshData_ext.G_BCs_VI_act      = G_BCs_VI_act;
    meshData_ext.G_BCs_VI_pas      = G_BCs_VI_pas;
    meshData_ext.G_CS_VI_act       = G_CS_VI_act;
    meshData_ext.G_CS_VI_pas       = G_CS_VI_pas;
    meshData_ext.G_flux_VI_pas     = G_flux_VI_pas;
    meshData_ext.G_B_pickup_VI_pas = G_B_pickup_VI_pas;
    meshData_ext.G_flux_VI_act     = G_flux_VI_act;
    meshData_ext.G_B_pickup_VI_act = G_B_pickup_VI_act;
    meshData_ext.vec_A             = vec_A;
    meshData_ext.n_cut_pas         = length(meshData_ext.ind_passive_cut);
    meshData_ext.n_G_source        = meshData_pas.n_G_source;
    meshData_ext.P_G_source        = meshData_pas.P_G_source;
    meshData_ext.w_G_source        = meshData_pas.w_G_source;
    
else
    
    fprintf('\n')
    disp('Computing VI matrices ...')
    
    [meshData_pas,...
        G_BCs_VI_act, ...
        G_BCs_VI_pas, ...
        G_flux_VI_pas, ...
        G_B_pickup_VI_pas, ...
        G_flux_VI_act, ...
        G_B_pickup_VI_act, ...
        vec_A] = fun_compute_VI_matrices_FRIDA_static(meshData_ext,meshData_pla,SETTINGS);
    
    n_pas = size(G_BCs_VI_pas,2);
    meshData_ext.n_pas = n_pas;
    
    meshData_ext.G_BCs_VI_act      = G_BCs_VI_act;
    meshData_ext.G_BCs_VI_pas      = G_BCs_VI_pas;
    meshData_ext.G_flux_VI_pas     = G_flux_VI_pas;
    meshData_ext.G_B_pickup_VI_pas = G_B_pickup_VI_pas;
    meshData_ext.G_flux_VI_act     = G_flux_VI_act;
    meshData_ext.G_B_pickup_VI_act = G_B_pickup_VI_act;
    meshData_ext.vec_A             = vec_A;
end


%% Save Data
fprintf('\n')
fprintf('\n')
fprintf('Saving data after preprocessing ... \n')

if SETTINGS.IS_EVOL
    save(['./data_in_FRIDA/INPUT_FRIDA_geo_preproc_', SETTINGS.filename, '.mat'], ...
        'meshData_pla', ...
        'meshData_ext',...
        'meshData_pas')
else
    save(['./data_in_FRIDA/INPUT_FRIDA_geo_preproc_', SETTINGS.filename, '.mat'], ...
        'meshData_pla', ...
        'meshData_ext')
end

pause(.01);
















