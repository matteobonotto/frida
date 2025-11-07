%% MATLAB 2 GMSH x RFX MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Write ITER geometrical data write in a suitable format for GMSH mesh
%   generation.
%   
%   Matteo Bonotto 01/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear;  clc;  close all;
restoredefaultpath

addpath ../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/
addpath ../../source/routines_FRIDA/


%% Input Geometry
input_Data = '';
load([input_Data ,'ITER_fw.mat']);
load([input_Data ,'ITER_geo_act.mat']);

%%
% % nt=size(meshData.t,1);
% % P1=[meshData.n(meshData.t(:,1),1) meshData.n(meshData.t(:,1),2)];
% % P2=[meshData.n(meshData.t(:,2),1) meshData.n(meshData.t(:,2),2)];
% % P3=[meshData.n(meshData.t(:,3),1) meshData.n(meshData.t(:,3),2)];
% % centro_t=(P1+P2+P3)/3;
% % figure
% %     hold on;    grid on
% % % %     pdemesh(p,[],t)
% %     mesh=triplot(meshData.t,meshData.n(:,1),meshData.n(:,2),'k');
% %     axis equal,  hold on, xlabel('r [m]'), ylabel('z [m]')
% %     axis([0.8 4.5 -1.5 1.5])
% % % Conduttori
% % ind_conduttori=find(keyreg>0);
% % %     triplot(meshData.t(ind_conduttori,:),meshData.n(:,1),meshData.n(:,2),'r')
% % % shell 
% % ind_shell=find(keyreg>=117 & keyreg<=176);
% %     triplot(meshData.t(ind_shell,:),meshData.n(:,1),meshData.n(:,2),'g')
% % % vessel
% % ind_vessel=find(keyreg>=57 & keyreg<=116);
% %     triplot(meshData.t(ind_vessel,:),meshData.n(:,1),meshData.n(:,2),'y')
% % % TSS
% % ind_TSS=find(keyreg>=177 & keyreg<=236);
% %     triplot(meshData.t(ind_TSS,:),meshData.n(:,1),meshData.n(:,2),'b')
% % % %     plot(centro_t(ind_TSS,1),centro_t(ind_TSS,2),'bo')
% % % saddle
% % ind_saddle=find(keyreg>=237 & keyreg<=244);
% %     triplot(meshData.t(ind_saddle,:),meshData.n(:,1),meshData.n(:,2),'c')
% % % %     plot(centro_t(ind_saddle,1),centro_t(ind_saddle,2),'co')
% % % Aactive coils (256 conduttori collegati in 194 circuiti)
% % ind_active=find(keyreg>=1 & keyreg<=56);
% %     triplot(meshData.t(ind_active,:),meshData.n(:,1),meshData.n(:,2),'r')
% %     
% %     
% %     
% % %%
% % 
% % % % for ii=237:244
% % % %     
% % % %    points = meshData.n(meshData.t(keyreg == ii,:),:);
% % % %    plot(points(:,1),points(:,2),'o')
% % % %    pause
% % % %     
% % % % end
% % for ii=237:244
% %     
% %    points = meshData.n(meshData.t(keyreg == ii,:),:);
% %    ind = boundary(points(:,1),points(:,2),0);
% %    points = points(ind(1:end-1),:);
% %    plot(points(:,1),points(:,2),'o', 'LineWidth',2)
% % % %    pause
% %     
% %    saddle_RFX{ii-236} = points;
% %    
% % end
% %     
% %     
% % % % save saddle_RFX saddle_RFX
% %     
% %     
% %     
% % %%
% % %     
% % %% VV, Shell,
% % 
% % rVV_int = 0.475;
% % rVV_est = 0.505;
% % 
% % rshell_int = 0.5115;
% % rshell_est = 0.5145;
% % 
% % rTSS_int = 0.553;
% % rTSS_est = 0.600;
% % 
% % nn_VV = 50;
% % nn_shell = 180;
% % nn_TSS = 80;
% % 
% % 
% % VV_int = [1.995+rVV_int*cos(linspace(0,2*pi-2*pi/nn_VV,nn_VV))' ...
% %     rVV_int*sin(linspace(0,2*pi-2*pi/nn_VV,nn_VV))'];
% % VV_ext = [1.995+rVV_est*cos(linspace(0,2*pi-2*pi/nn_VV,nn_VV))' ...
% %     rVV_est*sin(linspace(0,2*pi-2*pi/nn_VV,nn_VV))'];
% % 
% % shell_int = [1.995+rshell_int*cos(linspace(0,2*pi-2*pi/nn_shell,nn_shell))' ...
% %     rshell_int*sin(linspace(0,2*pi-2*pi/nn_shell,nn_shell))'];
% % shell_ext = [1.995+rshell_est*cos(linspace(0,2*pi-2*pi/nn_shell,nn_shell))' ...
% %     rshell_est*sin(linspace(0,2*pi-2*pi/nn_shell,nn_shell))'];
% % 
% % TSS_int = [1.995+rTSS_int*cos(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))' ...
% %     rTSS_int*sin(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))'];
% % TSS_ext = [1.995+rTSS_est*cos(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))' ...
% %     rTSS_est*sin(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))'];
% % 
% % 
% % figure
% % plot(VV_int(:,1),VV_int(:,2),'or-','LineWidth',2); hold on; axis equal;
% % plot(VV_ext(:,1),VV_ext(:,2),'or-','LineWidth',2)
% % % % plot(shell_int(:,1),shell_int(:,2),'ob-','LineWidth',2)
% % % % plot(shell_ext(:,1),shell_ext(:,2),'ob-','LineWidth',2)
% % plot(TSS_int(:,1),TSS_int(:,2),'og-','LineWidth',2)
% % plot(TSS_ext(:,1),TSS_ext(:,2),'og-','LineWidth',2)
% % 
% % delta_gap = pi/180/2;
% % theta_1 = linspace(0,pi-delta_gap,floor(nn_shell/2));
% % theta_2 = linspace(pi-delta_gap,2*pi,floor(nn_shell/2));
% % theta_2 = theta_2(1:end-1);
% % theta_3 = linspace(pi+delta_gap,pi-delta_gap+2*pi,floor(nn_shell)) - 2*pi;
% % 
% % VV_new = [1.995 0] + [rshell_est*cos(theta_3).' rshell_est*sin(theta_3).'; ...
% %     rshell_int*cos(fliplr(theta_3)).' rshell_int*sin(fliplr(theta_3)).'];
% % 
% % plot(VV_new(:,1),VV_new(:,2),'o-b','LineWidth',2)
% % 
% % % % figure
% % % % fill([VV_ext([1:end 1],1); VV_int([1:end 1],1)],[VV_ext([1:end 1],2); VV_int([1:end 1],2)],'r', 'EdgeColor', 'none')
% % % % hold on; axis equal
% % % % fill([shell_ext([1:end 1],1); shell_int([1:end 1],1)],[shell_ext([1:end 1],2); shell_int([1:end 1],2)],'r', 'EdgeColor', 'none')
% % % % fill([TSS_ext([1:end 1],1); TSS_int([1:end 1],1)],[TSS_ext([1:end 1],2); TSS_int([1:end 1],2)],'r', 'EdgeColor', 'none')
% % % % 
% % % % savefig('RFX_passive_fill')
% % 
% % % % Ne=20;
% % % % Ne = floor(.15*n_FW);
% % % % r0=1.995;
% % % % z0=0;
% % % % a=.12;
% % % % b=0.12;
% % % % theta=linspace(0,2*pi,Ne)';
% % % % ell.RR=r0+a*cos(theta(1:end-1));
% % % % ell.ZZ=z0+b*sin(theta(1:end-1));
% % % % ellisse = [ell.RR ell.ZZ];
% % % % 
% % % % [ellisse_r,ellisse_z]=equispaced(ellisse(:,1),ellisse(:,2),Ne);
% % % % ellisse = [ellisse_r ellisse_z];
% % % % 
% % % % plot(ellisse(:,1),ellisse(:,2),'-k')

%% Geometry: conductors
conduttori={};
ii=1;
for j=1:12
% %     pause
    conduttori{ii}=[source_act.RR(j,:)' source_act.ZZ(j,:)'] ;
    ii=ii+1;
end
n_act = length(conduttori);
% conduttori={} ii=1->4 saddle



%%


clear GEOMETRY
for ii=1:n_act
    GEOMETRY.shape{ii}=conduttori{ii};
end


for ii=1:n_act
    GEOMETRY.spacing{ii}=.25;
end


for ii=1:n_act
    GEOMETRY.surface{ii}=ii;
end

N_order = 1;
automesh_path = '../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/';
[p,e,t]=autoMESH_complex('tmp_ITER',GEOMETRY,N_order,automesh_path);



figure
triplot(t(:,1:3),p(:,1),p(:,2),'k')
axis equal; hold on;


meshData_all.t = t;
meshData_all.n = p;
meshData.t = t;
meshData.n = p;
meshData.nt = size(t,1);
meshData.nn = size(p,1);

meshData.N_order = N_order;


%% Define connection between circuits (via KONNAX matrix)

KONNAX_ACT = eye(n_act);
meshData.KONNAX_ACT = KONNAX_ACT;


%% Define element type for active conductors (1 -> n_act)

centro_t = [(p(t(:,1),1)+p(t(:,2),1)+p(t(:,3),1))/3 ...
    (p(t(:,1),2)+p(t(:,2),2)+p(t(:,3),2))/3];

type = zeros(meshData.nt,1);

figure
hold on; axis equal

for ii = 1:n_act
    
    col = rand(3,1);
    conduttori_ii = conduttori{ii};
    ind = find(inpolygon(centro_t(:,1),centro_t(:,2),conduttori_ii(:,1),conduttori_ii(:,2)));
    triplot(meshData.t(ind,1:3),meshData.n(:,1),meshData.n(:,2), 'color', col)
    type(ind) = ii;
    
end


ind_act = 1:max(type);
meshData.ind_act = ind_act;

meshData.type = type;

meshData.ind_sources = ind_act;
meshData.ind_pas = [];
meshData.ind_passive_cut = [];

qq = numel(unique([meshData.t(meshData.type == 101,:); ...
    meshData.t(meshData.type == 102,:); ...
    meshData.t(meshData.type == 103,:)]));
fprintf('Number of passive nodes = %i\n', qq)
fprintf('Number of total   nodes = %i\n', meshData.nn)


%%


figure
hold on; axis equal
for ii = 1:n_act
   
    ind_el = find(KONNAX_ACT(ii,:));
    ind_tri = find(sum(meshData.type.' == ind_el.'));
    triplot(meshData.t(ind_tri,1:3),meshData.n(:,1),meshData.n(:,2),'color',rand(3,1));
% %     pause
    
end



%% Saving Data

meshData_ext = meshData;
if N_order == 1
    save(sprintf('%sITERlike_mesh_VI_linear.mat',input_Data), 'meshData_ext')
% %     save RFXmod_noTSS_asCerma0nl_mesh_VI_linear meshData
% % elseif N_order == 2
% %     save RFXmod_mesh_VI_quadratic meshData
% % elseif N_order == 3
% %     save RFXmod_mesh_VI_cubic meshData
end





%%
clear;  clc;  close all;
restoredefaultpath

% % addpath ../../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/
dir_FRIDA = '../../source/routines_FRIDA/';
addpath(genpath(dir_FRIDA))

input_Data = '';


%% geometry per evol (VI)

%%% Plasma domain (PDE solution)
load([input_Data ,'ITER_fw.mat']);
load([input_Data ,'ITER_VV.mat']);

shell_int = shell_out(1:2:end,:);   

figure
plot(shell_int(:,1),shell_int(:,2),'or-')
axis equal; hold on;

n_FW = 300;  
n_bc = ceil(.75*n_FW);

FW = fw;

bc=[];

n_FW = n_FW;
n_bc = n_bc;
factor_bc = .2;
N_order = 2;
MEX_OPT = true;

spacing = zeros(2,1);
spacing(1) = .3;
spacing(2) = .3;

dir_AutoMESH = [dir_FRIDA, '/AutoMESH_ver_2/'];
tic
[meshData_pla] = fun_buildmesh_pla_FRIDA_evol(...
    FW,n_FW,bc,n_bc,factor_bc,spacing,N_order,dir_AutoMESH,MEX_OPT);
toc

%% Define a coupling surface (here we are taking the BCs curve as CS)
tri = meshData_pla.t;
ind_n_B = meshData_pla.ind_n_bc;
[ind_t_B] = fun_indTriOnEdge(tri,ind_n_B);
ind_surf = zeros(size(ind_t_B));

for ii = 1 : length(ind_t_B)
    tri_ii = meshData_pla.t(ind_t_B(ii),:);
    tmp = intersect(tri_ii,meshData_pla.ind_n_bc);
    ind_surf(ii) = tmp(end);
end

C_surf = meshData_pla.n(ind_surf,:);
meshData_pla.C_surf = C_surf;

plot(meshData_pla.C_surf(:,1), meshData_pla.C_surf(:,2), '*b')
legend('vac + plasma mesh', ...
    'plasma mesh', ...
    'limiter curve', ...
    'BCs curve', ...
    'normal unit vec',...
    'coupling surface')



%% Save mesh
save([input_Data, 'ITERlike_mesh_pla.mat'], 'meshData_pla')


return






