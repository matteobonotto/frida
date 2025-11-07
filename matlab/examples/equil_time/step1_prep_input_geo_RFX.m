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

dir_FRIDA = '../../source/routines_FRIDA/';
dir_AutoMESH = '../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/';
addpath(genpath(dir_FRIDA))

%% Input Geometry
input_Data = '../../data/geo/';
load([input_Data ,'RFX_create_l_saddle_01.mat'],'keyreg');
load([input_Data ,'RFX_create_l_saddle_01_linear.mat']);
load([input_Data ,'kon.mat']);

%%
nt=size(meshData.t,1);
P1=[meshData.n(meshData.t(:,1),1) meshData.n(meshData.t(:,1),2)];
P2=[meshData.n(meshData.t(:,2),1) meshData.n(meshData.t(:,2),2)];
P3=[meshData.n(meshData.t(:,3),1) meshData.n(meshData.t(:,3),2)];
centro_t=(P1+P2+P3)/3;
figure
    hold on;    grid on
% %     pdemesh(p,[],t)
    mesh=triplot(meshData.t,meshData.n(:,1),meshData.n(:,2),'k');
    axis equal,  hold on, xlabel('r [m]'), ylabel('z [m]')
    axis([0.8 4.5 -1.5 1.5])
% Conduttori
ind_conduttori=find(keyreg>0);
%     triplot(meshData.t(ind_conduttori,:),meshData.n(:,1),meshData.n(:,2),'r')
% shell 
ind_shell=find(keyreg>=117 & keyreg<=176);
    triplot(meshData.t(ind_shell,:),meshData.n(:,1),meshData.n(:,2),'g')
% vessel
ind_vessel=find(keyreg>=57 & keyreg<=116);
    triplot(meshData.t(ind_vessel,:),meshData.n(:,1),meshData.n(:,2),'y')
% TSS
ind_TSS=find(keyreg>=177 & keyreg<=236);
    triplot(meshData.t(ind_TSS,:),meshData.n(:,1),meshData.n(:,2),'b')
% %     plot(centro_t(ind_TSS,1),centro_t(ind_TSS,2),'bo')
% saddle
ind_saddle=find(keyreg>=237 & keyreg<=244);
    triplot(meshData.t(ind_saddle,:),meshData.n(:,1),meshData.n(:,2),'c')
% %     plot(centro_t(ind_saddle,1),centro_t(ind_saddle,2),'co')
% Aactive coils (256 conduttori collegati in 194 circuiti)
ind_active=find(keyreg>=1 & keyreg<=56);
    triplot(meshData.t(ind_active,:),meshData.n(:,1),meshData.n(:,2),'r')
    
    
    
%%

% % for ii=237:244
% %     
% %    points = meshData.n(meshData.t(keyreg == ii,:),:);
% %    plot(points(:,1),points(:,2),'o')
% %    pause
% %     
% % end
for ii=237:244
    
   points = meshData.n(meshData.t(keyreg == ii,:),:);
   ind = boundary(points(:,1),points(:,2),0);
   points = points(ind(1:end-1),:);
   plot(points(:,1),points(:,2),'o', 'LineWidth',2)
% %    pause
    
   saddle_RFX{ii-236} = points;
   
end
    
    
% % save saddle_RFX saddle_RFX
    
    
    
%%
%     
%% VV, Shell,

rVV_int = 0.475;
rVV_est = 0.505;

rshell_int = 0.5115;
rshell_est = 0.5145;

rTSS_int = 0.553;
rTSS_est = 0.600;

nn_VV = 50;
nn_shell = 180;
nn_TSS = 80;


VV_int = [1.995+rVV_int*cos(linspace(0,2*pi-2*pi/nn_VV,nn_VV))' ...
    rVV_int*sin(linspace(0,2*pi-2*pi/nn_VV,nn_VV))'];
VV_ext = [1.995+rVV_est*cos(linspace(0,2*pi-2*pi/nn_VV,nn_VV))' ...
    rVV_est*sin(linspace(0,2*pi-2*pi/nn_VV,nn_VV))'];

shell_int = [1.995+rshell_int*cos(linspace(0,2*pi-2*pi/nn_shell,nn_shell))' ...
    rshell_int*sin(linspace(0,2*pi-2*pi/nn_shell,nn_shell))'];
shell_ext = [1.995+rshell_est*cos(linspace(0,2*pi-2*pi/nn_shell,nn_shell))' ...
    rshell_est*sin(linspace(0,2*pi-2*pi/nn_shell,nn_shell))'];

TSS_int = [1.995+rTSS_int*cos(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))' ...
    rTSS_int*sin(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))'];
TSS_ext = [1.995+rTSS_est*cos(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))' ...
    rTSS_est*sin(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))'];


figure
plot(VV_int(:,1),VV_int(:,2),'or-','LineWidth',2); hold on; axis equal;
plot(VV_ext(:,1),VV_ext(:,2),'or-','LineWidth',2)
% % plot(shell_int(:,1),shell_int(:,2),'ob-','LineWidth',2)
% % plot(shell_ext(:,1),shell_ext(:,2),'ob-','LineWidth',2)
plot(TSS_int(:,1),TSS_int(:,2),'og-','LineWidth',2)
plot(TSS_ext(:,1),TSS_ext(:,2),'og-','LineWidth',2)

delta_gap = pi/180/2;
theta_1 = linspace(0,pi-delta_gap,floor(nn_shell/2));
theta_2 = linspace(pi-delta_gap,2*pi,floor(nn_shell/2));
theta_2 = theta_2(1:end-1);
theta_3 = linspace(pi+delta_gap,pi-delta_gap+2*pi,floor(nn_shell)) - 2*pi;

VV_new = [1.995 0] + [rshell_est*cos(theta_3).' rshell_est*sin(theta_3).'; ...
    rshell_int*cos(fliplr(theta_3)).' rshell_int*sin(fliplr(theta_3)).'];

plot(VV_new(:,1),VV_new(:,2),'o-b','LineWidth',2)

% % figure
% % fill([VV_ext([1:end 1],1); VV_int([1:end 1],1)],[VV_ext([1:end 1],2); VV_int([1:end 1],2)],'r', 'EdgeColor', 'none')
% % hold on; axis equal
% % fill([shell_ext([1:end 1],1); shell_int([1:end 1],1)],[shell_ext([1:end 1],2); shell_int([1:end 1],2)],'r', 'EdgeColor', 'none')
% % fill([TSS_ext([1:end 1],1); TSS_int([1:end 1],1)],[TSS_ext([1:end 1],2); TSS_int([1:end 1],2)],'r', 'EdgeColor', 'none')
% % 
% % savefig('RFX_passive_fill')

% % Ne=20;
% % Ne = floor(.15*n_FW);
% % r0=1.995;
% % z0=0;
% % a=.12;
% % b=0.12;
% % theta=linspace(0,2*pi,Ne)';
% % ell.RR=r0+a*cos(theta(1:end-1));
% % ell.ZZ=z0+b*sin(theta(1:end-1));
% % ellisse = [ell.RR ell.ZZ];
% % 
% % [ellisse_r,ellisse_z]=equispaced(ellisse(:,1),ellisse(:,2),Ne);
% % ellisse = [ellisse_r ellisse_z];
% % 
% % plot(ellisse(:,1),ellisse(:,2),'-k')

%% Geometria Saddle coils
conduttori={};
% % for ii=1:numel(ind_saddle)
% %     triplot(meshData.t(ind_saddle(ii),:),meshData.n(:,1),meshData.n(:,2),'k')
% %     pause
% % end
ii=1;
for j=1:4
    ind_saddle_coil=ind_saddle((j-1)*4+1:j*4);
    ind_n_saddle=reshape(meshData.t(ind_saddle_coil,:),...
        numel(meshData.t(ind_saddle_coil,:)),1);
    ind_n_saddle=unique(ind_n_saddle);
    coils.rr=meshData.n(ind_n_saddle,1);
    coils.zz=meshData.n(ind_n_saddle,2);
    [coils_sort]=fun_ordinapunti([coils.rr coils.zz]);
    plot(coils.rr,coils.zz,'.r')
% %     pause
    conduttori_saddle{ii}=coils_sort;
    ii=ii+1;
end
% conduttori={} ii=1->4 saddle


ii=1;
for j=237:244
    ind_saddle_coil = find(keyreg == j);
    ind_n_saddle=reshape(meshData.t(ind_saddle_coil,:),...
        numel(meshData.t(ind_saddle_coil,:)),1);
    ind_n_saddle=unique(ind_n_saddle);
    coils.rr=meshData.n(ind_n_saddle,1);
    coils.zz=meshData.n(ind_n_saddle,2);
    [coils_sort]=fun_ordinapunti([coils.rr coils.zz]);
    plot(coils_sort([1:end 1],1),coils_sort([1:end 1],2),'-or', 'LineWidth',2)
% %     pause
    conduttori_saddle_correct{ii}=coils_sort;
    ii=ii+1;
end


%% Geometria Attivi
load([input_Data ,'RFX_create_l_saddle_01.mat']); 
meshData.t=t(1:6,:)';
meshData.n=p';
ii=5;


load([input_Data ,'RFX_Geometria_Conduttori.mat']);
for ii=1:size(cond,2)
    conduttore=cond{ii};
    plot(conduttore(:,1),conduttore(:,2),'-*k')
% %     pause
    conduttori_attivi{ii}=cond{ii};
end


% % load([input_Data ,'RFX_Geometria_Conduttori.mat']);
% % for ii=1:size(cond,2)
% %     conduttore=cond{ii};
% %     
% %     if ii == 51+4
% %        conduttore = [2.386 .541; ...
% %            2.44 .541; ...
% %            2.44 .669; 
% %            2.386 .669];
% %     elseif ii == 52+4
% %        conduttore = [2.386 -.541; ...
% %            2.44 -.541; ...
% %            2.44 -.669; 
% %            2.386 -.669];
% %     elseif ii == 53+4
% %        conduttore = [2.567 0.343; ...
% %            2.621 0.343; ...
% %            2.621 0.471; ...
% %            2.567 0.471];    
% %     elseif ii == 54+4
% %        conduttore = [2.567 -0.343; ...
% %            2.621 -0.343; ...
% %            2.621 -0.471; ...
% %            2.567 -0.471];    
% %     end
% %     hold on
% % % %         plot3(conduttore(:,1),0*conduttore(:,2),conduttore(:,2),'-or', 'LineWidth',2)
% %     plot(conduttore(:,1),conduttore(:,2),'-or', 'LineWidth',2)
% %     
% % % %     pause
% %     conduttori_attivi{ii}=conduttore;
% % end


%%


clear GEOMETRY
GEOMETRY.shape{1} = VV_int;
GEOMETRY.shape{2} = VV_ext;
GEOMETRY.shape{3} = shell_int;
GEOMETRY.shape{4} = shell_ext;

GEOMETRY.shape{5} = TSS_int;

for ii=1:numel(conduttori_saddle)
    GEOMETRY.shape{ii+5}=conduttori_saddle{ii};
end
GEOMETRY.shape{10} = TSS_ext;
% % 
offset = numel(GEOMETRY.shape);
for ii=1:numel(conduttori_attivi)
    GEOMETRY.shape{ii+offset}=conduttori_attivi{ii};
end




for ii=1:10
    GEOMETRY.spacing{ii}=.1;
end


for ii=11:numel(GEOMETRY.shape)
    GEOMETRY.spacing{ii}=.05;
end

GEOMETRY.surface{1} = [2 1];
GEOMETRY.surface{2} = [4 3];
for ii=1:4
    GEOMETRY.surface{ii+2}=ii+5;
end
GEOMETRY.surface{7} = [10 5 6 7 8 9];
% % 
offset = numel(GEOMETRY.surface);
for ii=1:numel(conduttori_attivi)
    GEOMETRY.surface{ii+offset}=ii+10;
end

N_order = 1;
automesh_path = '../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/';
[p,e,t]=autoMESH_complex('tmp_RFX',GEOMETRY,N_order,automesh_path);



figure
triplot(t(:,1:3),p(:,1),p(:,2),'k')
axis equal; hold on;
plot(p(:,1),p(:,2),'*b', 'LineWidth',2)
plot(VV_int(:,1),VV_int(:,2),'or-','LineWidth',2); hold on; axis equal;
plot(VV_ext(:,1),VV_ext(:,2),'or-','LineWidth',2)
plot(shell_int(:,1),shell_int(:,2),'ob-','LineWidth',2)
plot(shell_ext(:,1),shell_ext(:,2),'ob-','LineWidth',2)
plot(TSS_int(:,1),TSS_int(:,2),'og-','LineWidth',2)
plot(TSS_ext(:,1),TSS_ext(:,2),'og-','LineWidth',2)


r_temp = 0.610;

temp = [1.995+r_temp*cos(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))' ...
    r_temp*sin(linspace(0,2*pi-2*pi/nn_TSS,nn_TSS))'];


numel(find(inpolygon(p(:,1),p(:,2),temp(:,1),temp(:,2))))


meshData_all.t = t;
meshData_all.n = p;




meshData.t = t;
meshData.n = p;
meshData.nt = size(t,1);
meshData.nn = size(p,1);

meshData.N_order = N_order;


%% KONNAX_ACT

% % type_attivi_OH = (1:5*8)';
% % type_attivi_FS = [41 41 42 42 43 43 44 44 45:56]';
% % 
% % type_attivi = [type_attivi_OH; type_attivi_FS];
% % 
% % type_saddle = [57 58 59 60]';
% % 
% % type_attivi_all = [type_attivi; type_saddle];
% % 
% % 
% % KONNAX_ACT = zeros(14,max(type_attivi_all));
% % 
% % % Ohmic Heating
% % type_attivi_OH_MAT = reshape(type_attivi_OH,5,8)';
% % % % temp_rows = [repmat((1:4)',1,5); repmat((1:4)',1,5)];
% % temp_rows = repmat([1 1 2 2 3 3 4 4]',1,5);
% % 
% % for ii=1:4
% %     KONNAX_ACT(temp_rows(ii,:),type_attivi_OH_MAT(ii,:)) = 1;
% %     KONNAX_ACT(temp_rows(ii+4,:),type_attivi_OH_MAT(ii+4,:)) = 1;
% % end
% % 
% % figure
% % spy(KONNAX_ACT)
% % 
% % % Field Shaping
% % jj_rows = 1;
% % jj_cols = 1;
% % for ii = max(type_attivi_OH)+1:2:max(type_attivi_FS)
% %     
% % % %     if ismember(ii,[41 42 43 44])
% % % %         KONNAX_ACT(max(temp_rows(:))+jj_rows,max(type_attivi_OH_MAT(:))+[jj_cols:jj_cols+3]) = 1;
% % % %         jj_rows = jj_rows+1;
% % % %         jj_cols = jj_cols+4;
% % % %         
% % % %     else
% %         KONNAX_ACT(max(temp_rows(:))+jj_rows,max(type_attivi_OH_MAT(:))+[jj_cols jj_cols+1]) = 1;
% %         jj_rows = jj_rows+1;
% %         jj_cols = jj_cols+2;
% %         
% % % %     end
% %     
% % end
% % 
% % spy(KONNAX_ACT)
% % 
% % 
% % % Saddle
% % KONNAX_ACT(13, [57 58]) = 1;
% % KONNAX_ACT(14, [59 60]) = -1;
% % % % KONNAX_ACT(13, [61 62]) = 1;
% % % % KONNAX_ACT(14, [63 64]) = -1;
% % spy(KONNAX_ACT)


type_attivi_OH = (1:5*8)';
type_attivi_FS = [41 41 42 42 43 43 44 44 45:56]';

type_attivi = [type_attivi_OH; type_attivi_FS];

type_saddle = [57:64]';

type_attivi_all = [type_attivi; type_saddle];


KONNAX_ACT = zeros(14,max(type_attivi_all));

% Ohmic Heating
type_attivi_OH_MAT = reshape(type_attivi_OH,5,8)';
% % temp_rows = [repmat((1:4)',1,5); repmat((1:4)',1,5)];
temp_rows = repmat([1 1 2 2 3 3 4 4]',1,5);

for ii=1:4
    KONNAX_ACT(temp_rows(ii,:),type_attivi_OH_MAT(ii,:)) = 1;
    KONNAX_ACT(temp_rows(ii+4,:),type_attivi_OH_MAT(ii+4,:)) = 1;
end

% % figure
% % spy(KONNAX_ACT)

% Field Shaping
jj_rows = 1;
jj_cols = 1;
for ii = max(type_attivi_OH)+1:2:max(type_attivi_FS)
    
    % %     temp = [find(type_attivi_FS == ii); find(type_attivi_FS == ii+1)];
    
    % %     KONNAX_ACT(max(temp_rows(:))+jj,max(type_attivi_OH_MAT(:))+temp) = 1;
    
    KONNAX_ACT(max(temp_rows(:))+jj_rows,max(type_attivi_OH_MAT(:))+[jj_cols jj_cols+1]) = 1;
    jj_rows = jj_rows+1;
    jj_cols = jj_cols+2;
    
end

% % spy(KONNAX_ACT)


% Saddle
KONNAX_ACT(13, [57 60 61 64]) = [1 1 -1 -1];
KONNAX_ACT(14, [58 59 62 63]) = [-1 1 1 -1];

meshData.KONNAX_ACT = KONNAX_ACT;


%%



centro_t = [(p(t(:,1),1)+p(t(:,2),1)+p(t(:,3),1))/3 ...
    (p(t(:,1),2)+p(t(:,2),2)+p(t(:,3),2))/3];


type = zeros(meshData.nt,1);

offset = 0;

figure
hold on; axis equal

for ii = 1:length(type_attivi_all)
    
    col = rand(3,1);
    if ii<=60
        conduttori_ii = conduttori_attivi{ii};
    else
            conduttori_ii = conduttori_saddle_correct{ii-60};
    end
    ind = find(inpolygon(centro_t(:,1),centro_t(:,2),conduttori_ii(:,1),conduttori_ii(:,2)));
    triplot(meshData.t(ind,1:3),meshData.n(:,1),meshData.n(:,2), 'color', col)
% %         pause
    type(ind) = type_attivi_all(ii);
    
end


% % figure; hold on; axis equal;
% % for ii = 1:max(type)
% % 
% %     int_t = find(type == ii);
% %     triplot(meshData.t(int_t,1:3), meshData.n(:,1), meshData.n(:,2), 'color' ,rand(1,3));
% %     pause
% %  
% % end

% % 
% % for ii = 41:2:48
% %     
% %     conduttori_ii = conduttori_attivi{ii};
% %     ind = find(inpolygon(centro_t(:,1),centro_t(:,2),conduttori_ii(:,1),conduttori_ii(:,2)));
% %     triplot(meshData.t(ind,1:3),meshData.n(:,1),meshData.n(:,2))
% %     % %     pause
% %     type(ind) = ii + offset;
% %     
% %     conduttori_ii = conduttori_attivi{ii+1};
% %     ind = find(inpolygon(centro_t(:,1),centro_t(:,2),conduttori_ii(:,1),conduttori_ii(:,2)));
% %     triplot(meshData.t(ind,1:3),meshData.n(:,1),meshData.n(:,2))
% %     % %     pause
% %     type(ind) = ii + offset;
% %     
% % end
% % 
% % for ii = 49:size(conduttori_attivi,2)
% %     
% %     conduttori_ii = conduttori_attivi{ii};
% %     ind = find(inpolygon(centro_t(:,1),centro_t(:,2),conduttori_ii(:,1),conduttori_ii(:,2)));
% %     triplot(meshData.t(ind,1:3),meshData.n(:,1),meshData.n(:,2))
% %     % %     pause
% %     type(ind) = ii + offset;
% % end
% % 
% % offset = size(conduttori_attivi,2);
% % 
% % for ii = 1:size(conduttori_saddle,2)
% %    
% %     conduttori_ii = conduttori_saddle{ii};
% %     ind = find(inpolygon(centro_t(:,1),centro_t(:,2),conduttori_ii(:,1),conduttori_ii(:,2)));
% %     triplot(meshData.t(ind,1:3),meshData.n(:,1),meshData.n(:,2),'r')
% % % %     pause
% %     type(ind) = ii + offset;
% % end

ind_act = 1:max(type);
meshData.ind_act = ind_act;



ind_VV = setdiff(find(inpolygon(centro_t(:,1),centro_t(:,2),VV_ext(:,1),VV_ext(:,2))), ...
    find(inpolygon(centro_t(:,1),centro_t(:,2),VV_int(:,1),VV_int(:,2))));
ind_shell = setdiff(find(inpolygon(centro_t(:,1),centro_t(:,2),shell_ext(:,1),shell_ext(:,2))), ...
    find(inpolygon(centro_t(:,1),centro_t(:,2),shell_int(:,1),shell_int(:,2))));
ind_TSS = setdiff(find(inpolygon(centro_t(:,1),centro_t(:,2),TSS_ext(:,1),TSS_ext(:,2))), ...
    find(inpolygon(centro_t(:,1),centro_t(:,2),TSS_int(:,1),TSS_int(:,2))));
ind_TSS = setdiff(ind_TSS,find(sum(type.' == (57:64).')));

type(ind_VV) = 101;
type(ind_shell) = 102;
type(ind_TSS) = 103;

triplot(meshData.t(ind_VV,1:3),meshData.n(:,1),meshData.n(:,2),'g')
triplot(meshData.t(ind_shell,1:3),meshData.n(:,1),meshData.n(:,2),'k')
triplot(meshData.t(ind_TSS,1:3),meshData.n(:,1),meshData.n(:,2),'c')



ind_pas = 101:103;
meshData.ind_pas = ind_pas;


figure
triplot(t(:,1:3),p(:,1),p(:,2),'k')
axis equal; hold on;
for ii = 1:max(type)
    try
        triplot(t(type == ii,1:3),p(:,1),p(:,2),'color',rand(3,1))
    catch
    end
    % %     pause
end

meshData.type = type;

ind_passive_cut = [2 3];

rho_VV    = 8.85e-6; % CREATE-L
rho_shell = 0.018e-6;
rho_TSS   = 1.08e-6;
meshData.rho_pas = [rho_VV; rho_shell; rho_TSS];


meshData.ind_sources = ind_act;
meshData.ind_pas = ind_pas;
meshData.ind_passive_cut = ind_passive_cut;

qq = numel(unique([meshData.t(meshData.type == 101,:); ...
    meshData.t(meshData.type == 102,:); ...
    meshData.t(meshData.type == 103,:)]));
fprintf('Number of passive nodes = %i\n', qq)
fprintf('Number of total   nodes = %i\n', meshData.nn)


% % figure
% % axis equal; hold on
% % triplot(meshData.t(find(sum(type.' == (101:103).')),1:3),meshData.n(:,1),meshData.n(:,2),'r')
% % triplot(meshData.t(find(sum(type.' == (1:60).')),1:3),meshData.n(:,1),meshData.n(:,2),'k')
% % box on;
% % 
% % savefig('RFX_passive_active_mesh')
% % 
% % figure
% % axis equal; hold on

%%


figure
hold on; axis equal
for ii = 1:14
   
    ind_el = find(KONNAX_ACT(ii,:));
    ind_tri = find(sum(meshData.type.' == ind_el.'));
    triplot(meshData.t(ind_tri,1:3),meshData.n(:,1),meshData.n(:,2),'color',rand(3,1));
% %     pause
    
end



%% Saving Data

meshData_ext = meshData;
if N_order == 1
    save(sprintf('%sRFXmod_mesh_VI_linear.mat',input_Data), 'meshData_ext')
% %     save RFXmod_noTSS_asCerma0nl_mesh_VI_linear meshData
% % elseif N_order == 2
% %     save RFXmod_mesh_VI_quadratic meshData
% % elseif N_order == 3
% %     save RFXmod_mesh_VI_cubic meshData
end

find(type == 57)




%%
clear;  clc;  close all;
restoredefaultpath


dir_FRIDA = '../../source/routines_FRIDA/';
dir_AutoMESH = '../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/';
addpath(genpath(dir_FRIDA))
input_Data = '../../data/geo/';


%%

%% geometry per evol (VI)

%%% Plasma domain (PDE solution)

%RFX-mod2
rfw = 0.49; 
r_PSS_int = 0.5115;
rbc = .5*(r_PSS_int + rfw);

%RFX-mod
rfw = 0.459; 
r_VV_int = 0.475;
rbc = .5*(r_VV_int + rfw);



nn_FW = 120;

FW = [1.995+rfw*cos(linspace(0,2*pi-2*pi/nn_FW,nn_FW))' ...
    rfw*sin(linspace(0,2*pi-2*pi/nn_FW,nn_FW))'];

bc=[1.995+rbc*cos(linspace(0,2*pi-2*pi/nn_FW,nn_FW))' ...
    rbc*sin(linspace(0,2*pi-2*pi/nn_FW,nn_FW))'];

n_FW = [];
n_bc = nn_FW;
factor_bc = [];
N_order = 2;
MEX_OPT = true;

spacing = zeros(2,1);
spacing(1) = .03;
spacing(2) = .03;

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
save([input_Data, 'RFXmod_mesh_pla.mat'], 'meshData_pla')


return






