clear
close all
clc
restoredefaultpath


%%


filename_save = '_ITER';
filename_preproc = '_ITER';


%% SETTINGS

% % addpath ../routines_FEMBEM_ver_1.2/
% % addpath ../routines_FRIDA_ver_1.9/
addpath ../routines_FRIDA_ver_2.1/

% % addpath ../routines_FEMBEM_ver_1/

% % filename_save = '_ITER_00090'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_HHR_00135'; filename_preproc = '_HHR_ITER';
% % filename_save = '_ITER_HR_00135'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_MR_00135'; filename_preproc = '_MR_ITER';


% % filename_save = '_ITER_MR_00083'; filename_preproc = '_MR_ITER';
filename_save = '_ITER_MR_00051'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_MR_00155'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_MR_00048'; filename_preproc = '_MR_ITER';
% % filename_save = '_FRIDA_evol_MR_00135'; filename_preproc = '_FRIDA_evol_MR_00135';

% % filename_save = '_ITER_LR_00090'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_MR_00083'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_MR_00155'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_LR_00051'; filename_preproc = '_LR_ITER';
% % filename_save = '_ITER_MR_00155'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_LR_00155'; filename_preproc = '_LR_ITER';
% % filename_save = '_ITER_MR_00135'; filename_preproc = '_MR_ITER';
% % filename_save = '_ITER_HR_00135'; filename_preproc = '_HR_ITER';

filename_out = filename_save;




%% SETTINGS

SETTINGS.filename_save = filename_save;
SETTINGS.filename_out = filename_out;
SETTINGS.filename_preproc = filename_preproc;

SETTINGS.PREPROC     = 'T';
SETTINGS.POSTPROC    = 'F';
SETTINGS.SAVE_OUTPUT = 'F';

SETTINGS.PREPROC_REORDERING = 'F';

SETTINGS.FIRST_SOLUTION = 'F';

SETTINGS.NN_SIZE_MAX_PARALLEL_JACOBIAN = 100;

SETTINGS.N_GAUSS = 12;

SETTINGS.shape_functions_N_order = 2;

% % SETTINGS.VAC_METHOD=1; % Vacuum solution with full Integral formulation
SETTINGS.VAC_METHOD=2; % Vacuum solution with FEM and integral BC

SETTINGS.FIGURES= 1;           % 0 -> no figures while iterating, 1 -> display figs
SETTINGS.FIGURES_DEBUG= 0;     % 0 -> no figures while iterating, 1 -> display figs

% % SETTINGS.plaparameter.Centroid = [2.05 -0.05]; % ok con 39136
% % SETTINGS.plaparameter.Centroid = [2.02 -0.05]; % ok con 39068
% % SETTINGS.plaparameter.Centroid = [3.1 0]; % 

% % SETTINGS.INITIAL_PLASMA_MODEL = 1; % Single filament 
SETTINGS.INITIAL_PLASMA_MODEL = 2; % filamentary model
% % SETTINGS.INITIAL_PLASMA_MODEL = 3; % jphi with Gaussian distribution

% % SETTINGS.SOLVER='PICARD';
SETTINGS.SOLVER='NR';
% % SETTINGS.SOLVER='NK';

SETTINGS.NEWTONRRAPHSON.RELAXED = 'F';

SETTINGS.NEWTONRRAPHSON.HMSP = 'T';
SETTINGS.NEWTONRRAPHSON.HMSP_THRESHOLD = .1;


SETTINGS.NEWTONRRAPHSON.TRESHOLD_1 = 15;
SETTINGS.NEWTONRRAPHSON.TRESHOLD_2 = 50;
SETTINGS.NEWTONRRAPHSON.TRESHOLD_3 = 150;
SETTINGS.NEWTONRRAPHSON.TRESHOLD_4 = 500;

SETTINGS.NEWTONRRAPHSON.ALPHA_1 = 1;
SETTINGS.NEWTONRRAPHSON.ALPHA_2 = .5;
SETTINGS.NEWTONRRAPHSON.ALPHA_3 = .2;
SETTINGS.NEWTONRRAPHSON.ALPHA_4 = .05;
SETTINGS.NEWTONRRAPHSON.ALPHA_5 = .005;


SETTINGS.QuasiNR = 'F';
SETTINGS.QuasiNR_fullsteps = 3;
SETTINGS.QuasiNR_step = 3;

% % % % SETTINGS.JACOBIAN=1;  % (standard, more reliable) semi-analytical computation of Jacobian
% % SETTINGS.JACOBIAN=2; % (fast, less reliable) semi-analytical computation of Jacobian
% % % % SETTINGS.JACOBIAN=3; % full-analytical computation of Jacobian

SETTINGS.SAVE_OUTPUT_FIG=0;
SETTINGS.SAVE_OUTPUT_MAT=0;

SETTINGS.FULLEQUIL=0;

SETTINGS.TOLL   = 5e-20;
SETTINGS.SCARTO = 1e-15; %(del residuo)

SETTINGS.ii_freeB_max=20;

SETTINGS.THRESHOLD_FIXB = 5E-3;
SETTINGS.TOLL_FIXB      = 2E-7; %2E-7;
SETTINGS.N_ITERMAX_FIXB = 50; %50

tolleranza_fix = SETTINGS.TOLL_FIXB;
N_ITERMAX_FIXB  = SETTINGS.N_ITERMAX_FIXB;
SETTINGS.NFREEB_FIXB = 10;

SETTINGS.perturtation = sqrt(eps);
% % SETTINGS.perturtation = 1e-6;

SETTINGS.J_PARAMETRIZATION.TYPE = 1; % use profiles of f(psi) and p(psi)
% % SETTINGS.J_PARAMETRIZATION.TYPE = 2; % use Blum's parametrization

SETTINGS.J_PARAMETRIZATION.INTERP1_METHOD = 'linear';
% % SETTINGS.J_PARAMETRIZATION.INTERP1_METHOD = 'pchip';

SETTINGS.RUN_MEX_ROUTINE = 'T';

% Smoothing profiles
SETTINGS.SMOOTH_FDFPROFILE = 'F';
SETTINGS.SMOOTH_FDFPROFILENPT = 7;

SETTINGS.SMOOTH_DPPROFILE = 'F';
SETTINGS.SMOOTH_DPPROFILENPT = 7;


SETTINGS.POSTPROCESSING_BSENS = 'F';
SETTINGS.POSTPROCESSING_PSISENS = 'F';
SETTINGS.FIELD_NODES = 'F';
SETTINGS.FIELD_GAUSS = 'F';
SETTINGS.POSTPROCESSING_BGRID = 'F';


% Postprocessing
SETTINGS.POSTPROCESSING_NGRIDR = 200;
SETTINGS.POSTPROCESSING_NGRIDZ = 200;


%% Run code
[OUTPUT_FEMBEM]=run_FRIDA_static(SETTINGS);


%% Post Processing

solk1             = OUTPUT_FEMBEM.solk1;
residuo_tot       = OUTPUT_FEMBEM.residuo_tot;
residuo_norm_tot  = OUTPUT_FEMBEM.residuo_norm_tot;
scarto_tot        = OUTPUT_FEMBEM.scarto_tot;
meshData          = OUTPUT_FEMBEM.meshData;
meshData_loc      = OUTPUT_FEMBEM.meshData_loc;
index             = OUTPUT_FEMBEM.index;
index_loc         = OUTPUT_FEMBEM.index_loc;
plaparameter      = OUTPUT_FEMBEM.plaparameter;
POSTPROC          = OUTPUT_FEMBEM.POSTPROC;



Rplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,1))...
    max(meshData_loc.n(meshData_loc.ind_n_bc,1))]+[-0.3 0.3];
Zplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,2))...
    max(meshData_loc.n(meshData_loc.ind_n_bc,2))]+[-0.3 0.3];




rr=meshData_loc.n(:,1);
zz=meshData_loc.n(:,2);
rgrid=linspace(Rplot(1),Rplot(2),200);
zgrid=linspace(Zplot(1),Zplot(2),200);
[RR,ZZ] = meshgrid(rgrid,zgrid);


% % Rplot=[1.4 2.6];
% % Zplot=[-.6 .6];

gg_pla = solk1.Jphi;



tic
F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),solk1.Jphi);
J_interp = NaN(meshData_loc.nn,1);
J_interp(index_loc.ind_n_INpla) = F(meshData_loc.n(index_loc.ind_n_INpla,1),meshData_loc.n(index_loc.ind_n_INpla,2));
toc


Jpla.faces=meshData_loc.t(:,1:3);
Jpla.vertices=meshData_loc.n;
Jpla.facevertexcdata=J_interp;


fontsize = 11;


figure

subplot(1,3,1)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
contour(RR,ZZ,PSI,60);
colormap('jet'); 
contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'-r','LineWidth',2);
% % plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-r','LineWidth',2);
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
title('$\psi$','Interpreter','latex','fontsize',fontsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize)
xlabel('r','Interpreter','latex','fontsize',fontsize); 
ylabel('z','Interpreter','latex','fontsize',fontsize); 

subplot(1,3,2)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
% % contour(RR,ZZ,gg_pla_grid,60);
hh=patch(Jpla,'facecolor','interp','edgecolor','none');
axis equal;  hold on;% title('j_\phi');
axis([Rplot Zplot])
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-r','LineWidth',2);
title('$j(\psi)$','Interpreter','latex','fontsize',fontsize);
xlabel('r','Interpreter','latex','fontsize',fontsize); 
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize)
% %       colormap('hot');


subplot(1,3,3)
% %       semilogy(scarto_tot,'r.-');
semilogy(1:numel(residuo_norm_tot),residuo_norm_tot,'k.-','LineWidth',2);
hold on
xlim([1 numel(residuo_norm_tot)])
title('$||F_k||/||x_k||$','Interpreter','latex','FontSize',fontsize)
grid on; axis square; box on
yticks([ 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1])
xticks([1 5 10 15 20])
xlabel('\# iterations','Interpreter','latex','fontsize',fontsize); 
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize)




ind_OUTPUT_fig = floor(1000*rand);
figure(ind_OUTPUT_fig)
subplot(2,2,1)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
hh=patch(Jpla,'facecolor','interp','edgecolor','none');
axis equal;  hold on; colorbar vert; % title('j_\phi');
axis([Rplot Zplot])
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-r');

% %       colormap('hot');

subplot(2,2,2)
% % H_K.faces=meshData_loc.t(index_loc.ind_t_Invess,1:3);
H_K.faces=meshData_loc.t(:,1:3);
H_K.vertices=meshData_loc.n;
H_K.facevertexcdata=solk1.h_k(1:end-3);
hh=patch(H_K,'facecolor','interp','edgecolor','none');
colorbar vert
axis square;  hold on;  colorbar vert;
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1)
axis equal; axis([Rplot Zplot]),  hold on; box on
axis([Rplot Zplot])

subplot(2,2,3)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
contour(RR,ZZ,PSI,60);
colormap('jet'); colorbar vert;
contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'k');
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'.r');
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')

subplot(2,2,4)
% %       semilogy(scarto_tot,'r.-');
semilogy(1:numel(residuo_tot),residuo_tot,'b.-'); hold on
semilogy(1:numel(residuo_norm_tot),residuo_norm_tot,'k.-');
semilogy(1:numel(scarto_tot),scarto_tot,'r.-');
xlim([1 numel(residuo_norm_tot)])
legend('||R_k||','||R_k||/||x_k||', 'scarto_k'); grid on; axis square
% %       xlim([1 100])
% %       ylim([min(residuo_tot) max(residuo_tot)])

%

fontize = 10;


figure

subplot(1,3,1)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
contour(RR,ZZ,PSI,60);
colormap('jet'); 
contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'-r','LineWidth',2);
% % plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-r','LineWidth',2);
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
% % axis([1.35 2.65 -.65 .65])
title('$\psi$','Interpreter','latex','fontsize',fontize);


subplot(1,3,2)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
% % contour(RR,ZZ,gg_pla_grid,60);
hh=patch(Jpla,'facecolor','interp','edgecolor','none');
axis equal;  hold on;% title('j_\phi');
axis([Rplot Zplot])
% % axis([1.35 3.65 -.65 .65])
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-r');
title('$j(\psi)$','Interpreter','latex','fontsize',fontize);
% %       colormap('hot');


subplot(1,3,3)
% %       semilogy(scarto_tot,'r.-');
semilogy(1:numel(residuo_norm_tot),residuo_norm_tot,'k.-','LineWidth',2);
xlim([1 numel(residuo_norm_tot)])
title('$||F_k||/||x_k||$','Interpreter','latex','FontSize',fontize)
grid on; box on; axis square; 
yticks([ 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1])
xticks([1 5 10 15 20])


% % initial_guess = [solk1.Psi; solk1.Psi_a; solk1.Psi_b; solk1.lambda];
% % save initial_guess initial_guess

rr=meshData_loc.n(:,1);
zz=meshData_loc.n(:,2);
rgrid=linspace(Rplot(1),Rplot(2),200);
zgrid=linspace(Zplot(1),Zplot(2),200);
[RR,ZZ] = meshgrid(rgrid,zgrid);


% % Rplot=[1.4 2.6];
% % Zplot=[-.6 .6];

gg_pla = zeros(size(solk1.Jphi));
gg_pla(index_loc.ind_Gauss_INpla) = solk1.Jphi(index_loc.ind_Gauss_INpla);

gg_pla_grid = griddata(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),gg_pla,RR,ZZ);


ind_OUTPUT_fig = floor(1000*rand);
figure(ind_OUTPUT_fig)
subplot(2,2,1)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
contour(RR,ZZ,gg_pla_grid,60);
axis equal;  hold on; colorbar vert; % title('j_\phi');
axis([Rplot Zplot])
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-k');
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')

% %       colormap('hot');

subplot(2,2,2)
% % H_K.faces=meshData_loc.t(index_loc.ind_t_Invess,1:3);
H_K.faces=meshData_loc.t(:,1:3);
H_K.vertices=meshData_loc.n;
H_K.facevertexcdata=solk1.h_k(1:end-3);
hh=patch(H_K,'facecolor','interp','edgecolor','none');
colorbar vert
axis square;  hold on;  colorbar vert;
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1)
axis equal; axis([Rplot Zplot]),  hold on; box on
axis([Rplot Zplot])

subplot(2,2,3)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
contour(RR,ZZ,PSI,60);
colormap('jet'); colorbar vert;
contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'k');
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'.r');
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')

subplot(2,2,4)
% %       semilogy(scarto_tot,'r.-');
semilogy(1:numel(residuo_tot),residuo_tot,'b.-'); hold on
semilogy(1:numel(residuo_norm_tot),residuo_norm_tot,'k.-');
semilogy(1:numel(scarto_tot),scarto_tot,'r.-');
xlim([1 numel(residuo_norm_tot)])
legend('||R_k||','||R_k||/||x_k||', 'scarto_k'); grid on; axis square



fontize = 10;

F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),solk1.Jphi);
J_interp = F(meshData_loc.n(:,1),meshData_loc.n(:,2));

JP.faces=meshData_loc.t(index_loc.ind_t_INpla,1:3);
% % JP.faces=meshData_loc.t(index_loc.ind_t_Invess,1:3);
JP.vertices=meshData_loc.n;
% % JP.facevertexcdata=NaN(meshData_loc.nn,1);
% % JP.facevertexcdata(index_loc.ind_n_INpla)=J_interp(index_loc.ind_n_INpla);
JP.facevertexcdata=J_interp;

figure
subplot(1,3,2)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
% % contour(RR,ZZ,gg_pla_grid,60);
hh=patch(JP,'facecolor','interp','edgecolor','none');
axis equal;  hold on;% title('j_\phi');
axis([2.5 10  Zplot]),  hold on,
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'-r','LineWidth',2);
title('$j(\Psi)$','Interpreter','latex','fontsize',fontize);
% %       colormap('hot');

subplot(1,3,1)
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([2.5 10  Zplot]),  hold on,
PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
contour(RR,ZZ,PSI,60);
colormap('jet'); 
% % contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'k');
contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'-r','LineWidth',2);
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
title('$\Psi$','Interpreter','latex','fontsize',fontize);

subplot(1,3,3)
% %       semilogy(scarto_tot,'r.-');
semilogy(1:numel(residuo_norm_tot),residuo_norm_tot,'k.-','LineWidth',2);
xlim([1 numel(residuo_norm_tot)])
title('$||R_k||/||x_k||$','Interpreter','latex','FontSize',fontize)
grid on; axis square; box on
yticks([ 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1])
xticks([1 5 10 15 20])



% % contour(h.XData,h.YData,h.ZData, [h.LevelList h.LevelList],'--k','LineWidth',2);



subplot(1,3,1); 
set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontize);
xlabel('$r$ [m]','Interpreter','latex','fontsize',fontize);
ylabel('$z$ [m]','Interpreter','latex','fontsize',fontize);
subplot(1,3,2); 
set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontize);
xlabel('$r$ [m]','Interpreter','latex','fontsize',fontize);
subplot(1,3,3); 
set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontize);
xlabel('Iterations','Interpreter','latex','fontsize',fontize);


namefig = sprintf('fig_articolo%s',num2str(filename_out));
% % print(namefig,'-depsc', '-r500')
% % savefig(namefig)




if SETTINGS.SAVE_OUTPUT_FIG==1
    figure(ind_OUTPUT_fig)
    name_OUTPUT_fig = sprintf('fig_OUT_%i_1', str2double(plaparameter.equilname));
    savefig(name_OUTPUT_fig);
end

try
    namefig = sprintf('fig_%s_CREATE', (plaparameter.equilname));
    openfig(namefig);
    hold on;
    edges=pdemesh(meshData.n',meshData.e',[]);
    set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
    axis([Rplot Zplot]),  hold on,
    PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
    % %       contour(RR,ZZ,PSI,100);
    % %       contour(RR,ZZ,PSI,[min(solk1.Psi_B) max(solk1.Psi_B)],'k');
    plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'.k');
    plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
    
    CATCH
end

if SETTINGS.SAVE_OUTPUT_FIG==1
    name_OUTPUT_fig = sprintf('fig_OUT_%i_2', str2double(plaparameter.equilname));
    savefig(name_OUTPUT_fig);
end

return







figure
plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.h_k(1:end-3),'o')


openfig('usn_create_BOUNDARY');
hold on; 
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,
PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1.Psi,RR,ZZ);
% %       contour(RR,ZZ,PSI,100);
      contour(RR,ZZ,PSI,[-.11855 -.11855],'r');
plot(solk1.Separatrix(:,1),solk1.Separatrix(:,2),'.k');

figure
hold on; 
edges=pdemesh(meshData.n',meshData.e',[]);
set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
axis([Rplot Zplot]),  hold on,

ref = load([plaparameter.input_Equil '39100_equilibrium08.mat'], 'R_boundary', 'Z_boundary');
ref = load([plaparameter.input_Equil '39122_equilibrium085AT.mat'], 'R_boundary', 'Z_boundary');



plot(ref.R_boundary,ref.Z_boundary)

