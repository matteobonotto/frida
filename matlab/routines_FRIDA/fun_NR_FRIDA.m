function [solk1]=fun_NR_FRIDA(INPUT_solver)

%%
mu0=4*pi*1.e-7;

%%

solk         = INPUT_solver.solk;
meshData     = INPUT_solver.meshData;
meshData_loc = INPUT_solver.meshData_loc;
plaparameter = INPUT_solver.plaparameter;
ii_freeB     = INPUT_solver.ii_freeB;
SETTINGS     = INPUT_solver.SETTINGS;   

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    Centroid = plaparameter.Centroid;
    Ipla =     plaparameter.Ipla;
    psibar =   plaparameter.psibar;
    FdF =      plaparameter.FdF;
    dP =      plaparameter.dP;
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    R_0 =      plaparameter.R_0;
    beta_0 =   plaparameter.beta_0;
    alpha_M =  plaparameter.alpha_M;
    alpha_N =  plaparameter.alpha_N;
    Ipla =     plaparameter.Ipla;
end


%%
nn              = meshData_loc.nn;
nn_D            = length(meshData_loc.ind_D);
nn_B            = length(meshData_loc.ind_B);
nt              = meshData_loc.nt;
N_order         = meshData_loc.shape_functions_N_order;
ww_Gauss        = meshData_loc.ww_Gauss_pla;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
shape_functions = meshData_loc.shape_functions;
Iphi            = zeros(nn,1);
ind_D           = meshData_loc.ind_D;
ind_B           = meshData_loc.ind_B;
KK_nobc         = meshData_loc.KK_nobc;
KK_D            = KK_nobc(ind_D,ind_D);
KK_bc           = KK_nobc(ind_D,ind_B);
GG_BEM          = meshData_loc.GG_BEM;
nodes_Gauss_pla = meshData_loc.nodes_Gauss_pla;
ww_Gauss_pla    = meshData_loc.ww_Gauss_pla;
rr              = meshData_loc.n(:,1);
rr_Gauss        = meshData_loc.nodes_Gauss_pla(:,1);
ind_n_INpla     = meshData_loc.ind_n_INpla;
% % meshData_loc.KK_D_tilde = -(GG_BEM.psi/mu0)*KK_D;
% % meshData_loc.KK_bc_tilde = (speye(nn_B)-(GG_BEM.psi/mu0)*KK_bc);


%%
[~,ind] = max(solk.Psi(solk.ind_n_axis));
solk.Psi_axis_node = solk.Psi(solk.ind_n_axis(ind));


%% Fixed boundary step to converge on current density

if SETTINGS.RUN_FIXBOUNDARY == true && ...
        solk.scarto_norm >= SETTINGS.THRESHOLD_FIXB && ...
        ii_freeB <= SETTINGS.NFREEB_FIXB
    
    fprintf('Fixed-Boundary step \n')
    
    input_FixB.plaparameter = plaparameter;
    input_FixB.meshData_loc = meshData_loc;
    input_FixB.solk         = solk;
    input_FixB.SETTINGS     = SETTINGS;
    
    tic
    [solk_fix] = fun_FixB_FRIDA(input_FixB);
    toc
    
    solk.Psi = solk_fix.Psi;
    solk.Psi_axis = solk_fix.Psi_axis;
    solk.ind_n_axis = solk_fix.ind_n_axis;
    solk.delta_a_weights = solk_fix.delta_a_weights;
    
    solk.lambda = solk_fix.lambda;
    solk.lambda_star = solk_fix.lambda_star;
    
    [~,ind] = max(solk.Psi(solk.ind_n_axis));
    % % solk.ind_n_axis = solk.ind_n_axis(ind);
    solk.Psi_axis_node = solk.Psi(solk.ind_n_axis(ind));
    % % solk.delta_a_weights = 1;
    
    clear solk_fix
    
else
    
    fprintf('Skipping Fixed-Boundary step \n')

% %     solk.lambda_star = mu0*solk.lambda;
    
end

    Rplot=[.9*min(meshData_loc.n(:,1)) 1.1*max(meshData_loc.n(:,1))];
    Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
    rgrid=linspace(Rplot(1),Rplot(2),200);
    zgrid=linspace(Zplot(1),Zplot(2),200);
    [RR,ZZ] = meshgrid(rgrid,zgrid);
    PSI_2 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk.Psi,RR,ZZ);
    % %
    
    if SETTINGS.FIGURES_DEBUG == 1
        figure
        edges=pdemesh(meshData.n',meshData.e_interface',[]);
        set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
        axis([Rplot Zplot]),  hold on,
        contour(RR,ZZ,PSI_2,100);
        contour(RR,ZZ,PSI_2,[solk.Psi_Bpla solk.Psi_Bpla],'r','LineWidth',2);
        colormap('cool'); colorbar vert;
        plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
    end
    
lambda_k = solk.lambda;
lambda_star_k = solk.lambda_star;


%%


% % tic
% % [Source,solk] = fun_calcSourceFEM_adaptive(meshData_loc,solk,plaparameter,SETTINGS);
% % toc
% % 
% % Source_D = Source(ind_D);
% % 
% % Psi_axis=solk.Psi_axis;
% % Psi_Bpla=sum(solk.Psi_B)./numel(solk.Psi_B);
% % gg_tot = solk.gg_tot;
% % Psi_k = solk.Psi_k;
% % 
% % if SETTINGS.J_PARAMETRIZATION_TYPE == 2
% %     aa_r_tot =  solk.aa_r_tot;
% % end
% % 

%%
tri             = meshData_loc.t;
N_order         = meshData_loc.shape_functions_N_order;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
n_Gauss         = meshData_loc.n_Gauss;
shape_functions = meshData_loc.shape_functions;

Psi_k=solk.Psi;
Psi_axis=solk.Psi_axis_node;
Psi_Bpla=sum(solk.Psi_Bpla)./numel(solk.Psi_Bpla);
psi_bar=(Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux

if SETTINGS.RUN_MEX_ROUTINE == false
    [psi_bar] = fun_evaluateFluxGaussPoints_v2(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [psi_bar] = fun_evaluateFluxGaussPoints_v2_mex(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
end
solk.psi_bar = psi_bar;



if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    
% %     psi_bar(psi_bar>1) = NaN;
    psi_bar(psi_bar<0) = 0;

    fdfn=interp1(psibar,FdF,psi_bar,'linear');   
    dpn=interp1(psibar,dP,psi_bar,'linear');   
    gg_tot=(rr_Gauss.*dpn+fdfn./(mu0*rr_Gauss)); %computation of nodes plasma current density
    gg_tot=gg_tot/2/pi;
    
    qq = isnan(gg_tot);
    gg_tot(isnan(gg_tot)) = 0;
    
    gg_tot_in = gg_tot;
    gg_tot_in(meshData_loc.ind_G_out) = 0;

    nn = meshData_loc.nn;
    ind_t_selected = meshData_loc.ind_t_INpla;
    ww_Gauss = meshData_loc.ww_Gauss_pla;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
        [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
    
    % %     figure;  plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),Source,'o')
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    
    tmp_gg = ((1-psi_bar.^alpha_M).^alpha_N);
    tmp_gg(imag(tmp_gg) ~= 0) = abs(tmp_gg(imag(tmp_gg) ~= 0)); % occhio!!!
    aa_r_tot = (rr_Gauss*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss);
    gg_tot = tmp_gg.*aa_r_tot;

    gg_tot_in = gg_tot;
    gg_tot_in(meshData_loc.ind_G_out) = 0;
    
    nn = meshData_loc.nn;
    ind_t_selected = meshData_loc.ind_t_INpla;
    ww_Gauss = meshData_loc.ww_Gauss_pla;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
        [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot_in,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
    
end

Source_D = Source(ind_D);


%%
% %         F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),gg_tot);
% % 
% %         J_interp = F(meshData_loc.n(:,1),meshData_loc.n(:,2));
% % 
% %         Jpla.faces=meshData_loc.t(:,[1 4 2 5 3 6]);
% %         Jpla.vertices=meshData_loc.n;
% %         Jpla.facevertexcdata=J_interp;
% % 
% % 
% %         figure
% %         axis([Rplot Zplot]),  hold on,
% %         hh=patch(Jpla,'facecolor','interp','edgecolor','none');
% %         axis equal;  hold on; colorbar vert; % title('j_\phi');
% %         axis([Rplot Zplot])
% %         plot(solk.Separatrix(:,1),solk.Separatrix(:,2),'k')
% %         
% %         Rplot=[.9*min(meshData_loc.n(:,1)) 1.08*max(meshData_loc.n(:,1))];
% %         Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
% %         rgrid=linspace(Rplot(1),Rplot(2),200);
% %         zgrid=linspace(Zplot(1),Zplot(2),200);
% %         [RR,ZZ] = meshgrid(rgrid,zgrid);
% %         axis([Rplot Zplot]),  hold on,
% %         PSI = griddata(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),gg_tot,RR,ZZ);
% %         figure
% %         contour(RR,ZZ,PSI,50);
% %         axis equal
% % 
% %         
% %         figure;  plot3(meshData_loc.n(:,1), ...
% %             meshData_loc.n(:,2),Source,'o')


%% Sigma rows of stiffness matrix

sigma_a_tmp = zeros(1,nn);
sigma_b_tmp = zeros(1,nn);

sigma_a_tmp(solk.ind_n_axis) = solk.delta_a_weights;
sigma_b_tmp(solk.ind_n_XP)   = solk.delta_b_weights;

if SETTINGS.FIGURES_DEBUG == 1
    plot(meshData_loc.n(solk.ind_n_axis,1),meshData_loc.n(solk.ind_n_axis,2),'r*','LineWidth',2)
    plot(meshData_loc.n(solk.ind_n_XP,1),meshData_loc.n(solk.ind_n_XP,2),'b*','LineWidth',2)
end

sigma_a = [sigma_a_tmp(ind_D) sigma_a_tmp(ind_B)];
sigma_b = [sigma_b_tmp(ind_D) sigma_b_tmp(ind_B)];


%% Build stiffness matrix and RHS

GG_BEM_pla = zeros(nn_B,nn);
GG_BEM_pla(:,ind_D) = GG_BEM.psi;
GG_BEM_pla(:,meshData_loc.ind_n_NOTpla) = 0;
GG_BEM_pla = GG_BEM_pla(:,ind_D);

KK_D_tilde = -(GG_BEM_pla/mu0)*KK_D;
KK_bc_tilde = sparse(speye(nn_B)-(GG_BEM_pla/mu0)*KK_bc);


AA = [KK_D KK_bc sparse(nn_D,2) -Source_D; ...
     KK_D_tilde KK_bc_tilde sparse(nn_B,3); ...
    sigma_a -1 0 0; ...
    sigma_b 0 -1 0; ...
    sparse(1,nn) 0 0 sum(Source_D)]; 

RHS = [zeros(nn_D,1); solk.Psi_MS_bc; 0; 0; mu0*Ipla];

xx_k = [solk.Psi(ind_D); solk.Psi(ind_B); solk.Psi_axis; solk.Psi_B; solk.lambda_star];
% % xx_k = [solk.Psi(ind_D); solk.Psi(ind_B); solk.Psi_axis; max(solk.Psi_Bpla); solk.lambda_star];

ind_plot = floor(1000*rand(1));
if SETTINGS.FIGURES_DEBUG == 1
    figure(ind_plot);  hold on;
% %     plot(plaparameter.Psi_CREATE_loc_B,'k') ; hold on;
    plot(solk.Psi(ind_B),'r--','LineWidth',2)
end

%% Jacobian (numeric) 

perturbation = SETTINGS.perturbation;

tic
fprintf('Computing Jacobian matrix (semi-analytic) \n')
calc_Jacobian_FRIDA_semianalytic_fast
toc

tmp = DRpsibc_Dpsi(:);

tmp_pos = tmp(tmp>0);
tmp_mean_pos = sum(tmp_pos)/numel(tmp_pos);
% % tmp_pos2 = tmp_pos;
% % tmp_pos2(tmp_pos< .1*tmp_mean_pos)=0;


if SETTINGS.FAR_FIELD_SPARS == true
    
    threshold = SETTINGS.FAR_FIELD_SPARS_THRESHOLD;
    
    for ii=1:size(DRpsibc_Dpsi,1)
        DRpsibc_Dpsi(ii,abs(DRpsibc_Dpsi(ii,:))<threshold*tmp_mean_pos) = 0;
    end
    
end


% % calc_Jacobian_FRIDA_semianalytic_fast_v2

JAC = [DRpsi_Dpsi DRpsi_Dpsibc DRpsi_Dpsia DRpsi_Dpsib DRpsi_Dlambda; ...
    DRpsibc_Dpsi DRpsibc_Dpsibc DRpsibc_Dpsia DRpsibc_Dpsib DRpsibc_Dlambda; ...
    sigma_a DRpsia_Dpsia DRpsia_Dpsib DRpsia_Dlambda; ...
    sigma_b DRpsib_Dpsia DRpsib_Dpsib DRpsib_Dlambda; ...
    DRlambda_Dpsi DRlambda_Dpsibc DRlambda_Dpsia DRlambda_Dpsib DRlambda_Dlambda];




%% Solution

RR_psi_k = AA*xx_k-RHS;

fprintf('Computing solution \n')
time_solver_start = tic;

if SETTINGS.RUN_MEX_ROUTINE == false
    
    h_k = mldivide(JAC,RR_psi_k); % usually faster for purely FEM problems (purely sparse Matrix)
    
elseif SETTINGS.RUN_MEX_ROUTINE == true
    
    h_k = cs_lusol (JAC,RR_psi_k,1,[]); % usually faster for FEM-BEM problems
    % %     opts.ordering = 'amd';
    % %     h_k = umfpack (JAC, '\', RR_psi_k,[],opts); % usually faster than cs_lusol
    
end

toc(time_solver_start)

time_solver = toc(time_solver_start);

fprintf('\n')
fprintf('Residue psi(end) = %1.3e \n',h_k(end-3))
fprintf('Residue psi_a    = %1.3e \n',h_k(end-2))
fprintf('Residue psi_b    = %1.3e \n',h_k(end-1))
fprintf('Residue lambda   = %1.3e \n',h_k(end))
fprintf('\n')

xx_k1 = xx_k - h_k;

alfa = 1;

Ipla_check = xx_k1(end)*sum(Source_D)/mu0;
error_Ipla = 100*abs((Ipla_check - plaparameter.Ipla))/plaparameter.Ipla;

if error_Ipla > 1e-3
    fprintf('Preliminary check on lambda*Iphi:          error=%6.4f\n',error_Ipla)
else
    fprintf('Preliminary check on lambda*Iphi:          error=%6.3e\n',error_Ipla)
end

if SETTINGS.SOLVER_RELAXED == true

    treshold_1 = SETTINGS.SOLVER_RELAX_TRESHOLD_1;
    treshold_2 = SETTINGS.SOLVER_RELAX_TRESHOLD_2;
    treshold_3 = SETTINGS.SOLVER_RELAX_TRESHOLD_3;
    treshold_4 = SETTINGS.SOLVER_RELAX_TRESHOLD_4;

    alfa_1 = SETTINGS.SOLVER_RELAX_ALPHA_1;
    alfa_2 = SETTINGS.SOLVER_RELAX_ALPHA_2;
    alfa_3 = SETTINGS.SOLVER_RELAX_ALPHA_3;
    alfa_4 = SETTINGS.SOLVER_RELAX_ALPHA_4;
    alfa_5 = SETTINGS.SOLVER_RELAX_ALPHA_5;

    if error_Ipla < treshold_1
        alfa = alfa_1;
        
    elseif error_Ipla >= treshold_1 && error_Ipla < treshold_2
        alfa = alfa_2;
        
    elseif error_Ipla >= treshold_2 && error_Ipla < treshold_3
        alfa = alfa_3;
        
    elseif error_Ipla >= treshold_3 && error_Ipla < treshold_4
        alfa = alfa_4;
        
    elseif error_Ipla >= treshold_4
        alfa = alfa_5;
        
    end
        
% %     if error_Ipla < SETTINGS.NEWTONRRAPHSON.STABLE.TRESHOLD_1
% %         alfa = 1;
% %         
% %     elseif error_Ipla >= 20 && error_Ipla < 50
% %         alfa = .75;
% %         
% %     elseif error_Ipla >= 50 && error_Ipla < 150
% %         alfa = .2;
% %         
% %     elseif error_Ipla >= 150 && error_Ipla < 500
% %         alfa = .05;
% %         
% %     elseif error_Ipla >= 500
% %         alfa = .01;
% %         
% %     end
end

xx_k1 = xx_k - alfa*h_k;
    
% % xx_k1 = AA\(RHS);

Psi_tot = zeros(meshData_loc.nn,1);
Psi_tot(ind_D) = xx_k1(1:nn_D);
Psi_tot(ind_B) = xx_k1(nn_D+1:nn);

lambda_star = xx_k1(end);
lambda = lambda_star/mu0;
Psi_a = xx_k1(end-2);
Psi_b = xx_k1(end-1);

Ipla_check = lambda*sum(Source_D);

error_Ipla = 100*abs((Ipla_check - plaparameter.Ipla))/plaparameter.Ipla;

if alfa ~= 1
    fprintf('Check on lambda*Iphi after regularization: error=%6.4f\n',error_Ipla)
end

Iphi = Source*lambda;

Jphi = gg_tot*lambda;

% % tic
% % F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),Jphi);
% % J_interp = zeros(meshData_loc.nn,1);
% % J_interp(index_loc.ind_n_INpla) = F(meshData_loc.n(index_loc.ind_n_INpla,1),meshData_loc.n(index_loc.ind_n_INpla,2));
% % toc
% % 
% % 
% % JP.faces=meshData_loc.t(:,1:3);
% % JP.vertices=meshData_loc.n;
% % JP.facevertexcdata=J_interp;
% % 
% % figure
% % edges=pdemesh(meshData.n',meshData.e',[]);
% % set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
% % axis([Rplot Zplot]),  hold on,
% % hh=patch(JP,'facecolor','interp','edgecolor','none');
% % axis equal;  hold on;% title('j_\phi');
% % axis([Rplot Zplot])
% % axis([1.35 2.65 -.65 .65])
% % plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')



%%
if SETTINGS.FIGURES_DEBUG == 1
    figure(ind_plot)
    plot((Psi_tot(ind_B)),'b:','LineWidth',2)
end
      
if SETTINGS.FIGURES_DEBUG == 1
    figure
    edges=pdemesh(meshData.n',meshData.e_interface',[]);
    set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
    Rplot=[.9*min(meshData_loc.n(:,1)) 1.08*max(meshData_loc.n(:,1))];
    Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
    rgrid=linspace(Rplot(1),Rplot(2),200);
    zgrid=linspace(Zplot(1),Zplot(2),200);
    [RR,ZZ] = meshgrid(rgrid,zgrid);
    axis([Rplot Zplot]),  hold on,
    PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),Psi_tot,RR,ZZ);
    contour(RR,ZZ,PSI,50);
    colormap('cool'); colorbar vert;
    contour(RR,ZZ,PSI,[Psi_b Psi_b],'k');
    plot(meshData_loc.n(solk.ind_n_axis,1),meshData_loc.n(solk.ind_n_axis,2),'r*','LineWidth',2)
    plot(meshData_loc.n(solk.ind_n_XP,1),meshData_loc.n(solk.ind_n_XP,2),'b*','LineWidth',2)
    plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k.')
end


%% Output

solk1.Psi_MS_bc   = solk.Psi_MS_bc;
solk1.psi_bar     = psi_bar;
solk1.Ipla_check  = Ipla_check;
solk1.RR_psi_k    = RR_psi_k;
% % solk1.xx_k        = [Psi_tot; Psi_a; Psi_b; lambda];
solk1.xx_k          = [Psi_tot([ind_D; ind_B]); Psi_a; Psi_b; lambda_star];
solk1.h_k_psi     = (Psi_tot - solk.Psi);
solk1.h_k         = h_k;
solk1.Psi         = Psi_tot;
solk1.lambda      = lambda;
solk1.lambda_star = lambda_star;
solk1.error_Ipla  = error_Ipla;
solk1.Iphi        = Iphi;
solk1.Jphi        = Jphi;
solk1.Psi_a       = Psi_a;
solk1.Psi_b       = Psi_b;
solk1.time_solver = time_solver;




