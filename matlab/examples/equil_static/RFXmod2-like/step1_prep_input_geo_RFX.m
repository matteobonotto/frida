%% MATLAB RFXmod2 Electromagnetic structures poloidal cross section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   D Abate October 2023
% April 2024: corretto raggio dei tondini r=3mm (prima era erroneamente
% r=6mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;  clc;  close all;

% % addpath(genpath('../../../source/'))

%% Input Geometry
dir_FRIDA = '../../../source/routines_FRIDA/';
addpath(genpath(dir_FRIDA))

dir_autoMESH = [dir_FRIDA, '/AutoMESH_VI_ver_1.0/'];

%% Set folders

dir_equil = '../../../data/equil/';
dir_sensosr = '../../../data/mod2_sensors_positions/';
dir_geo = '../../../data/geo/';
dir_in_FRIDA = './data_in_FRIDA/';

load([dir_geo ,'RFX_create_l_saddle_01.mat'],'keyreg');
load([dir_geo ,'RFX_create_l_saddle_01_linear.mat']);
load([dir_geo ,'kon.mat']);


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
    

%% RFXMOD2%  
R0=1.995;
Z0 = 0;
%
% R_vv=R0;
% r_outer_vv=0.505;
theta = linspace(0,2*pi,100);
theta = theta(1:end-1);

r_TSS_int= 0.555;		%interno TSS
r_TSS_ext = 0.600;		%raggio esterno TSS
r_TSS_ave = .5*(r_TSS_int+r_TSS_ext);

RVTSS = 2;
VTSSin(1,:) = RVTSS+r_TSS_int*cos(theta);
VTSSin(2,:) = r_TSS_int*sin(theta);
VTSSext(1,:) = RVTSS+r_TSS_ext*cos(theta);
VTSSext(2,:) = r_TSS_ext*sin(theta);

r_inner_shell=0.5115;
r_outer_shell=0.5145;
r_shell_ave=(r_inner_shell+r_outer_shell)/2; %raggio medio scocca
SHELL_in(1,:)=R0 + r_inner_shell*cos(theta);
SHELL_in(2,:)=Z0 + r_inner_shell*sin(theta);
SHELL_ext(1,:)=R0 + r_outer_shell*cos(theta);
SHELL_ext(2,:)=Z0 + r_outer_shell*sin(theta);

SHELL_in = SHELL_in';
SHELL_ext=SHELL_ext';
VTSSin=VTSSin';
VTSSext = VTSSext';

rfw = 0.4905;
FW(1,:)=R0+rfw*cos(theta);
FW(2,:)=Z0 +rfw*sin(theta);
FW = FW';

figure
plot(SHELL_in(:,1),SHELL_in(:,2),'r-','LineWidth',2); hold on; axis equal;
plot(SHELL_ext(:,1),SHELL_ext(:,2),'r-','LineWidth',2)
% % plot(shell_int(:,1),shell_int(:,2),'ob-','LineWidth',2)
% % plot(shell_ext(:,1),shell_ext(:,2),'ob-','LineWidth',2)
plot(VTSSin(:,1),VTSSin(:,2),'g-','LineWidth',2)
plot(VTSSext(:,1),VTSSext(:,2),'g-','LineWidth',2)

%gabbia: 6 tondini toroidali
theta_cage = linspace(0,2*pi,10);
theta_cage = theta_cage(1:end-1);
r_cage = 3e-3;
% theta_cage(1) = (180-23.57)*pi/180;
% theta_cage(2) = (180-27.86)*pi/180;
% theta_cage(3) = (180-36.43)*pi/180;
% theta_cage(4) = (180+23.57)*pi/180;
% theta_cage(5) = (180+27.86)*pi/180;
% theta_cage(6) = (180+36.43)*pi/180;
% r_theta_cage = 535e-3;

CAGE{1}(:,1)=1504.6e-3 + r_cage*cos(theta_cage');
CAGE{1}(:,2)=213.9e-3 + r_cage*sin(theta_cage');

CAGE{2}(:,1)=1522e-3 + r_cage*cos(theta_cage');
CAGE{2}(:,2)=250e-3 + r_cage*sin(theta_cage');

CAGE{3}(:,1)=1564.5e-3 + r_cage*cos(theta_cage');
CAGE{3}(:,2)=317.7e-3 + r_cage*sin(theta_cage');

CAGE{4}(:,1)=1504.6e-3 + r_cage*cos(theta_cage');
CAGE{4}(:,2)=-213.9e-3 + r_cage*sin(theta_cage');

CAGE{5}(:,1)=1522e-3 + r_cage*cos(theta_cage');
CAGE{5}(:,2)=-250e-3 + r_cage*sin(theta_cage');

CAGE{6}(:,1)=1564.5e-3 + r_cage*cos(theta_cage');
CAGE{6}(:,2)=-317.7e-3 + r_cage*sin(theta_cage');
[CAGE{1}(:,1) CAGE{1}(:,2)]
plot(CAGE{1}(:,1),CAGE{1}(:,2),'m', 'LineWidth',2)
plot(CAGE{2}(:,1),CAGE{2}(:,2),'m', 'LineWidth',2)
plot(CAGE{3}(:,1),CAGE{3}(:,2),'m', 'LineWidth',2)
plot(CAGE{4}(:,1),CAGE{4}(:,2),'m', 'LineWidth',2)
plot(CAGE{5}(:,1),CAGE{5}(:,2),'m', 'LineWidth',2)
plot(CAGE{6}(:,1),CAGE{6}(:,2),'m', 'LineWidth',2)

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
   % plot(coils.rr,coils.zz,'.r')
% %     pause
    conduttori_saddle{ii}=coils_sort;
    ii=ii+1;
end
% % conduttori={} ii=1->4 saddle


ii=1;
for j=237:244
    ind_saddle_coil = find(keyreg == j);
    ind_n_saddle=reshape(meshData.t(ind_saddle_coil,:),...
        numel(meshData.t(ind_saddle_coil,:)),1);
    ind_n_saddle=unique(ind_n_saddle);
    coils.rr=meshData.n(ind_n_saddle,1);
    coils.zz=meshData.n(ind_n_saddle,2);
    [coils_sort]=fun_ordinapunti([coils.rr coils.zz]);
%     plot(coils_sort([1:end 1],1),coils_sort([1:end 1],2),'-or', 'LineWidth',2)
% %     pause
    conduttori_saddle_correct{ii}=coils_sort;
    plot(conduttori_saddle_correct{ii}([1:end 1],1),conduttori_saddle_correct{ii}([1:end 1],2),'-b', 'LineWidth',2)
    ii=ii+1;
end


%% Geometria Attivi
load([dir_geo ,'RFX_create_l_saddle_01.mat']); 
meshData.t=t(1:6,:)';
meshData.n=p';

load([dir_geo ,'RFX_Geometria_Conduttori.mat']);
for ii=1:size(cond,2)
    conduttore=cond{ii};
    plot(conduttore(:,1),conduttore(:,2),'-k', 'LineWidth',2)
% %     pause
    conduttori_attivi{ii}=cond{ii};
end
grid on
xlabel('R (m)')
ylabel('Z (m)')
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)

save RFX_mod2_geometry SHELL_ext SHELL_in VTSSext VTSSin CAGE conduttori_attivi conduttori_saddle_correct

%% mesh
% % addpath ./AutoMESH_VI_ver_1.0/


clear GEOMETRY
GEOMETRY.shape{1} = SHELL_in;
GEOMETRY.shape{2} = SHELL_ext;
GEOMETRY.shape{3} = VTSSin;
for ii=1:numel(conduttori_saddle)
    GEOMETRY.shape{ii+3}=conduttori_saddle{ii};
end
GEOMETRY.shape{ii+3+1} = VTSSext; %8
for kk=1:6 %9
GEOMETRY.shape{ii+3+1+kk} = CAGE{kk}; 
end
% % 
offset = numel(GEOMETRY.shape);
for ii=1:numel(conduttori_attivi)
    GEOMETRY.shape{ii+offset}=conduttori_attivi{ii};
end

for ii=1:8 % PSS VTSS SC
    GEOMETRY.spacing{ii}=.05;
end
for ii=9:14 % cage
    GEOMETRY.spacing{ii}=.025;
end

for ii=(offset+1):numel(GEOMETRY.shape)
    GEOMETRY.spacing{ii}=.05;
end

% SC
% for ii=57:64
%     GEOMETRY.spacing{ii}=.01;
% end

GEOMETRY.surface{1} = [2 1]; %PSS
for ii=1:4
    GEOMETRY.surface{ii+1}=ii+3;
end
GEOMETRY.surface{6} = [8 3 4 5 6 7]; %VTSS (8,3) SC (4,5,6,7)
for ii=1:6 %CAGE
    GEOMETRY.surface{ii+6}=ii+8; %9,10,11,12,13,14
end
% % 
offset = numel(GEOMETRY.surface);
for ii=1:numel(GEOMETRY.shape)
    GEOMETRY.surface{ii+offset}=ii+14;
end

N_order = 1;
[p,e,t]=autoMESH_complex(...
    'tmp_RFX',GEOMETRY,N_order,dir_autoMESH);


nodes_matrix = [p(t(:,1),:) ...
    p(t(:,2),:) ...
    p(t(:,3),:)];

P1 = nodes_matrix(:,[1 2]);
P2 = nodes_matrix(:,[3 4]);
P3 = nodes_matrix(:,[5 6]);

a = P2-P1;
b = P3-P1;
res = .5*abs(a(:,1).*b(:,2) - a(:,2).*b(:,1));
[a,b] = min(res)


figure
triplot(t(:,1:3),p(:,1),p(:,2),'k')
axis equal; hold on;
triplot(t(b,1:3),p(:,1),p(:,2),'r', 'linewidth', 10)

%plot(p(:,1),p(:,2),'b', 'LineWidth',2)
plot(SHELL_in(:,1),SHELL_in(:,2),'r-','LineWidth',2); hold on; axis equal;
plot(SHELL_ext(:,1),SHELL_ext(:,2),'r-','LineWidth',2)
plot(VTSSext(:,1),VTSSext(:,2),'b-','LineWidth',2)
plot(VTSSin(:,1),VTSSin(:,2),'b-','LineWidth',2)
plot(CAGE{1}(:,1),CAGE{1}(:,2),'m', 'LineWidth',2)
plot(CAGE{2}(:,1),CAGE{2}(:,2),'m', 'LineWidth',2)
plot(CAGE{3}(:,1),CAGE{3}(:,2),'m', 'LineWidth',2)
plot(CAGE{4}(:,1),CAGE{4}(:,2),'m', 'LineWidth',2)
plot(CAGE{5}(:,1),CAGE{5}(:,2),'m', 'LineWidth',2)
plot(CAGE{6}(:,1),CAGE{6}(:,2),'m', 'LineWidth',2)


r_temp = 0.610;
nn_TSS=100
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


figure; hold on; axis equal;
for ii = 57:max(type)

    int_t = find(type == ii);
    triplot(meshData.t(int_t,1:3), meshData.n(:,1), meshData.n(:,2), 'color' ,rand(1,3));
    hold on
    %pause

end

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

% % figure
% % triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2),'b')
% % axis equal


% ind_CAGE = setdiff(find(inpolygon(centro_t(:,1),centro_t(:,2),CAGE{1}(:,1),CAGE{1}(:,2))), ...
%     find(inpolygon(centro_t(:,1),centro_t(:,2),VTSSin(:,1),VTSSin(:,2))));
for ii=1:length(CAGE)
ind_CAGE{ii} = find(inpolygon(centro_t(:,1),centro_t(:,2),CAGE{ii}(:,1),CAGE{ii}(:,2)));
end
ind_shell = setdiff(find(inpolygon(centro_t(:,1),centro_t(:,2),SHELL_ext(:,1),SHELL_ext(:,2))), ...
    find(inpolygon(centro_t(:,1),centro_t(:,2),SHELL_in(:,1),SHELL_in(:,2))));
ind_VTSS = setdiff(find(inpolygon(centro_t(:,1),centro_t(:,2),VTSSext(:,1),VTSSext(:,2))), ...
    find(inpolygon(centro_t(:,1),centro_t(:,2),VTSSin(:,1),VTSSin(:,2))));
ind_VTSS = setdiff(ind_VTSS,find(sum(type.' == (57:64).')));

for jj=1:size(ind_CAGE,2)
    type(ind_CAGE{jj}) = 95+jj;
end
type(ind_shell) = 102;
type(ind_VTSS) = 103;

triplot(meshData.t(ind_VTSS,1:3),meshData.n(:,1),meshData.n(:,2),'b')
triplot(meshData.t(ind_shell,1:3),meshData.n(:,1),meshData.n(:,2),'r')
%triplot(meshData.t(ind_CAGE,1:3),meshData.n(:,1),meshData.n(:,2),'m')


ind_pas = 96:103;
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

ind_passive_cut = [7 8]; %1:6 gabbia, 7 pss, 8 vtss

rho_CAGE    = 1.29e-6; % CREATE-L
rho_shell = 0.018e-6;
rho_VTSS   = 7.2e-7;%1.08e-6;
meshData.rho_pas = [rho_CAGE*ones(6,1); rho_shell; rho_VTSS];


meshData.ind_sources = ind_act;
meshData.ind_pas = ind_pas;
meshData.ind_passive_cut = ind_passive_cut;

figure
triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2),'k')
axis equal
hold on
triplot(meshData.t(meshData.type==ind_pas(1),1:3),meshData.n(:,1),meshData.n(:,2),'b')
triplot(meshData.t(meshData.type==ind_pas(2),1:3),meshData.n(:,1),meshData.n(:,2),'r')
triplot(meshData.t(meshData.type==ind_pas(3),1:3),meshData.n(:,1),meshData.n(:,2),'c')




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



%% sensors

sens = readtable([dir_sensosr, 'b3_arrays.csv'], 'NumHeaderLines', 1); % Load sensors
flux_loops = readtable([dir_sensosr, 'Tor_Vloop.xlsx'], 'NumHeaderLines', 1); % Load sensors
% sensors geometrical data
sens_mod2_0501 = table2array(sens(1:3:42, 3:end)); %% primo array toroidale
sens_mod2_1102 = table2array(sens(43:3:84, 3:end)); %% secondo array toroidale staggerato
Points_0501 = [sens_mod2_0501(:,1) sens_mod2_0501(:,2)];
Points_1102 = [sens_mod2_1102(:,1) sens_mod2_1102(:,2)];
theta_0501 = transpose(sens_mod2_0501(:,4));
theta_1102 = transpose(sens_mod2_1102(:,4));
thetarad_0501 = theta_0501*pi/180;
thetarad_1102 = theta_1102*pi/180;

%%%%%%%%%%% gabbia
theta_cage(1) = (180-23.57)*pi/180;
theta_cage(2) = (180-27.86)*pi/180;
theta_cage(3) = (180-36.43)*pi/180;
theta_cage(4) = (180+23.57)*pi/180;
theta_cage(5) = (180+27.86)*pi/180;
theta_cage(6) = (180+36.43)*pi/180;

clear CAGE
CAGE(1,1)=1504.6e-3;
CAGE(1,2)=213.9e-3;
CAGE(2,1)=1522e-3;
CAGE(2,2)=250e-3;
CAGE(3,1)=1564.5e-3;
CAGE(3,2)=317.7e-3;
CAGE(4,1)=1504.6e-3;
CAGE(4,2)=-213.9e-3;
CAGE(5,1)=1522e-3;
CAGE(5,2)=-250e-3;
CAGE(6,1)=1564.5e-3;
CAGE(6,2)=-317.7e-3;
%%%%%%%%%%%

% thetarad = [thetarad_0501,thetarad_1102,theta_cage];
% BSENS_R = [Points_0501(:,1);Points_1102(:,1);CAGE(:,1)];
% BSENS_Z = [Points_0501(:,2);Points_1102(:,2);CAGE(:,2)];
%BSENS_T = [-sin(thetarad)' cos(thetarad)'];
%pickup = [BSENS_R BSENS_Z BSENS_T];
thetarad = [thetarad_0501,thetarad_1102,theta_cage];

ngrid=20;
[RR_grid,ZZ_grid] = meshgrid(linspace(1.4,2.6,ngrid),linspace(-.55,.55,ngrid));
RR_grid = reshape(RR_grid,ngrid*ngrid,1);
ZZ_grid = reshape(ZZ_grid,ngrid*ngrid,1);

%1:28 pickup, 29:34 gabbia, 35: 134 regione plasma - direzione R
BSENS_R = [Points_0501(:,1);Points_1102(:,1);CAGE(:,1);RR_grid];
BSENS_Z = [Points_0501(:,2);Points_1102(:,2);CAGE(:,2);ZZ_grid];
DIR_R = ones(length(BSENS_R),2).*[1 0];
DIR_Z = ones(length(BSENS_Z),2).*[0 1];
pickup = [BSENS_R BSENS_Z DIR_R; BSENS_R BSENS_Z DIR_Z];
sensors.pickup = pickup;

fluxloops_mod2 = table2array(flux_loops(1:8, 3:end));
Points_fluxloops = [fluxloops_mod2(:,2) fluxloops_mod2(:,4)];
theta_fluxloops = transpose(fluxloops_mod2(:,3));

% RGM sensori di flusso
% R, Z coordinates of conductor section center
%       SFP_xA           R           Z
%       SFP_1A       1.239       0.035
%       SFP_1A       1.239       0.211
%       SFP_2A       1.334       0.343
%       SFP_2A       1.334       0.593
%       SFP_3A       1.377       0.645
%       SFP_3A       1.484       0.799
%       SFP_4A       1.631       1.001
%       SFP_4A       1.631       1.247
%
%       SFP_xB           R           Z
%       SFP_1B       1.239      -0.035
%       SFP_1B       1.239      -0.211
%       SFP_2B       1.334      -0.343
%       SFP_2B       1.334      -0.593
%       SFP_3B       1.377      -0.645
%       SFP_3B       1.484      -0.799
%       SFP_4B       1.631      -1.001
%       SFP_4B       1.631      -1.247

RGM=load([dir_sensosr, 'RGMPolFluxSensorCoord.txt']);
R_fluxloops = [Points_fluxloops(:,1);CAGE(:,1);RGM(:,1)];
Z_fluxloops = [Points_fluxloops(:,2);CAGE(:,2);RGM(:,2)];
flux_loops = [R_fluxloops Z_fluxloops];
sensors.flux_loops = flux_loops;


%% Saving Data
% % dir_geo = './data_in_FRIDA/';
% % mkdir(dir_geo)
% % if N_order == 1
% %     save([dir_geo, 'RFXmod2_mesh_VI_linear_cage_v2.mat'], 'meshData')
% % % %     save RFXmod_noTSS_asCerma0nl_mesh_VI_linear meshData
% % % % elseif N_order == 2
% % % %     save RFXmod_mesh_VI_quadratic meshData
% % % % elseif N_order == 3
% % % %     save RFXmod_mesh_VI_cubic meshData
% % end

%%
% % clear;  clc;  close all;
dir_FRIDA = '../../../source/routines_FRIDA/';
dir_autoMESH = [dir_FRIDA, '/AutoMESH_VI_ver_1.0/'];

% % restoredefaultpath
% % 
% % % % addpath ../../../source/routines_FRIDA/AutoMESH_VI_ver_1.0/
% % dir_FRIDA = '../../source/routines_FRIDA/';
% % addpath(genpath(dir_FRIDA))
% % 
% % dir_geo = '';

%% Geometry I: Plasma domain
%RFX-mod2
rfw = 0.4905;
r_PSS_int = 0.5115;
rbc = .5*(r_PSS_int + rfw);

% %RFX-mod
% rfw = 0.459; 
% r_VV_int = 0.475;
% rbc = .5*(r_VV_int + rfw);

nn_FW = 120;

FW = [1.995+rfw*cos(linspace(0,2*pi-2*pi/nn_FW,nn_FW))' ...
    rfw*sin(linspace(0,2*pi-2*pi/nn_FW,nn_FW))'];

bc=[1.995+rbc*cos(linspace(0,2*pi-2*pi/nn_FW,nn_FW))' ...
    rbc*sin(linspace(0,2*pi-2*pi/nn_FW,nn_FW))'];

n_FW = [];
n_bc = nn_FW;
factor_bc = [];
N_order = 2;
MEX_OPT = false;

Conductors_geo = [];
Conductors_index = [];

spacing = zeros(size(Conductors_geo,2)+2,1);
spacing(1) = .03;
spacing(2) = .03;

tic
[meshData_pla] = fun_buildmesh_pla_FRIDA_evol(...
    FW,n_FW,bc,n_bc,factor_bc,spacing,N_order,dir_autoMESH,MEX_OPT);
toc

%%% Coupling surface
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
if ~exist(dir_in_FRIDA,'dir')
    mkdir(dir_in_FRIDA)
end
save(...
    [dir_in_FRIDA, 'tmp_INPUT_FRIDA_geo.mat'], ...
    'meshData', ...
    'meshData_pla', ...
    'sensors')




