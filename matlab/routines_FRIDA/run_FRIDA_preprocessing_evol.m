

%

N_order = meshData_loc.shape_functions_N_order;

tri   = meshData_loc.t;
nodes = meshData_loc.n;
nn    = meshData_loc.nn;


%%

if SETTINGS.FIGURES_DEBUG == 1
    figure
    hold on;    grid on
    mesh_vac=triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2));
    set(mesh_vac,'color',.7*[1 1 1]);
    
    for ii = 1:Conductors.Nconductors
        mesh_vac=triplot(meshData.t(meshData.type == ii,1:3),meshData.n(:,1),meshData.n(:,2));
        set(mesh_vac,'color','r');
    end
    
    edges=pdemesh(meshData.n',meshData.e_interface',[]);
    set(edges,'color','k','LineWidth',1);  axis equal;
    axis equal;  hold on;  box on; xlabel('r [m]'), ylabel('z [m]')
end


%% Shape functions (for Integral and FEM)
tic

disp('Computing shape functions ...')

if SETTINGS.RUN_MEX_ROUTINE == false
    [shape_f,shape_f_norm] = fun_shape_functions_stable(tri,nodes,N_order);    
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [shape_f,shape_f_norm] = fun_shape_functions_stable_mex(tri,nodes,N_order);
end

meshData_loc.shape_functions = shape_f;
meshData_loc.shape_functions_norm = shape_f_norm;

toc



%% Gauss nodes for plasma and stiffness matrix computation
disp('Plasma domain: computing Gauss points for plasma (loc)...')
tic

% % n_Gauss = 3*N_order;
GAUSS_QUAD_DEGREE_PLA = SETTINGS.GAUSS_QUAD_DEGREE_PLA;

tri = meshData_loc.t;

nodes_Pla = [meshData_loc.n(tri(:,1),:) ...
    meshData_loc.n(tri(:,2),:) ...
    meshData_loc.n(tri(:,3),:)];

P1 = nodes_Pla(:,[1 2]);
P2 = nodes_Pla(:,[3 4]);
P3 = nodes_Pla(:,[5 6]);

if SETTINGS.RUN_MEX_ROUTINE == false
    [ww_Gauss_pla,nodes_Gauss_pla,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,GAUSS_QUAD_DEGREE_PLA);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [ww_Gauss_pla,nodes_Gauss_pla,n_Gauss] = fun_Gauss_points_triangle_stable_mex(P1,P2,P3,GAUSS_QUAD_DEGREE_PLA);    
end

meshData_loc.n_Gauss = n_Gauss;
meshData_loc.nodes_Gauss_pla = nodes_Gauss_pla;
meshData_loc.ww_Gauss_pla = ww_Gauss_pla;


%% Incidence Matrix between nodes and triangles
disp('Building incidence Matrix between nodes and triangles (loc)...')

tic

if SETTINGS.RUN_MEX_ROUTINE == false
    [MatInd_nodes_tri] = fun_MatInd_nodes_tri(N_order,tri,nn);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [MatInd_nodes_tri] = fun_MatInd_nodes_tri_fast_mex(N_order,tri,nn);
end

meshData_loc.MatInd_nodes_tri=MatInd_nodes_tri;

toc

%% Incidence Matrix between triangles and Gauss points
disp('Building incidence Matrix between triangles and Gauss points (loc)...')
tic

index_cols = 1:meshData_loc.nt*meshData_loc.n_Gauss;
index_rows = reshape(repmat(1:meshData_loc.nt,meshData_loc.n_Gauss,1),meshData_loc.nt*meshData_loc.n_Gauss,1);

MatInd_tri_Gauss = reshape(index_cols,meshData_loc.n_Gauss,meshData_loc.nt)';

meshData_loc.MatInd_tri_Gauss=MatInd_tri_Gauss;

toc


%% Incidence Matrix between nodes
disp('Building Incidence Matrix between nodes (loc)...')
tic

[MatInd_nodes] = fun_MatInd_nodes_v2(tri,nodes,MatInd_nodes_tri);
meshData_loc.MatInd_nodes=MatInd_nodes;

toc



%% Stiffness Matrix (FEM)
disp('Computing stiffness matrix (loc)...')
tic

meshData_loc.GAUSS_QUAD_DEGREE_STIFFMAT = SETTINGS.GAUSS_QUAD_DEGREE_STIFFMAT;

% Already with MEX options!!!
[KK_nobc] = fun_StiffMat_FEM(meshData_loc,SETTINGS.RUN_MEX_ROUTINE);

meshData_loc.KK_nobc = KK_nobc;
toc

%% Boundary-Stiffness matrix (BEM) 
disp('Computing Green Matrix for BEM (loc)...')

tic

clear source point
point.RR=meshData_loc.n(meshData_loc.ind_B,1);
point.ZZ=meshData_loc.n(meshData_loc.ind_B,2);
source.R=meshData_loc.n(meshData_loc.ind_D,1);
source.Z=meshData_loc.n(meshData_loc.ind_D,2);

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
        
        psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
        psi = psi/2/pi;
        
        GG_BEM.psi(:,ii) = psi;
        
    end
end

meshData_loc.GG_BEM = GG_BEM;
toc





%% Plasma response equation (coupling surface)

C_surf = meshData.C_surf;
n_pla_eq = size(C_surf,1);

centro = sum(C_surf)/size(C_surf,1);
R0 = centro(1);
Z0 = centro(2);

[theta_C,rho_C] = cart2pol(C_surf(:,1)-R0,C_surf(:,2)-Z0);

% define a target surface where the flux by I_p has to match the reference one
r_target = 1.01*max(rho_C);
C_target = [R0+r_target*cos(theta_C) r_target*sin(theta_C)];


clear source point
point.RR=C_target(:,1);
point.ZZ=C_target(:,2);
source.R=C_surf(:,1);
source.Z=C_surf(:,2);
source.current=ones(numel(source.R),1);
GG_eq_target  = fun_GreenMat_Loop_flux( source, point );
% %     GG_BEM_eq = fun_Field_Loop_flux(source,point);

clear source point
point.RR=C_target(:,1);
point.ZZ=C_target(:,2);
source.R=meshData_loc.n(:,1);
source.Z=meshData_loc.n(:,2);
source.current=ones(numel(source.R),1);
GG_pla_target = fun_GreenMat_Loop_flux(source,point);

T_Ipla_eq = GG_eq_target.psi\GG_pla_target.psi;

meshData.T_Ipla_eq = T_Ipla_eq;

T_Ipla_eq_D = zeros(n_pla_eq,nn);
T_Ipla_eq_D(:,meshData_loc.ind_D) = T_Ipla_eq(:,meshData_loc.ind_D); 
T_Ipla_eq_D = T_Ipla_eq_D(:,meshData_loc.ind_D);

meshData.T_Ipla_eq_D = T_Ipla_eq_D;



%% Save Data
fprintf('Saving data after preprocessing ... \n')

save(['INPUT_FRIDA_geo_preproc_' SETTINGS.filename_geo], ...
    'meshData_loc', ...
    'meshData')


pause(.1);
















