function [solk_fix] = fun_FixB_FRIDA(input_FixB)



%%
plaparameter = input_FixB.plaparameter;
meshData_loc = input_FixB.meshData_loc;
solk         = input_FixB.solk;
SETTINGS     = input_FixB.SETTINGS;

%%

ind_n_ONpla = meshData_loc.ind_n_ONpla;
ind_B = meshData_loc.ind_B;
ind_D = meshData_loc.ind_D;
ind_n_bc = meshData_loc.ind_n_bc;

KK_nobc = meshData_loc.KK_nobc;

tri             = meshData_loc.t;
N_order         = meshData_loc.shape_functions_N_order;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
n_Gauss         = meshData_loc.n_Gauss;
shape_functions = meshData_loc.shape_functions;
ww_Gauss = meshData_loc.ww_Gauss_pla;
nn = meshData_loc.nn;

% % ind_t_on_curve = meshData_loc.ind_t_on_curve;

%%

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


mu0=4*pi*1.e-7;

%%
% % ind_n_NOTpla = meshData_locind_n_NOTpla;

xx_k = solk.xx_k;

% % Psi_fix_old = xx_k(1:meshData_loc.nn);

Psi_fix_old = zeros(meshData_loc.nn,1);
Psi_fix_old(meshData_loc.ind_D) = xx_k(1:meshData_loc.nn_D);
Psi_fix_old(meshData_loc.ind_B) = xx_k(meshData_loc.nn_D+1:meshData_loc.nn);

solk_fix = solk;

ind_D = ind_D;
ind_B = ind_B;

ind_B_fix = [ind_n_ONpla(:); ind_B(:)];

Psi_B_fix = [Psi_fix_old(ind_n_ONpla); solk.Psi(ind_n_bc)];

[ind_B_fix,bb,~] = unique(ind_B_fix,'stable');
Psi_B_fix = Psi_B_fix(bb);

ind_D_fix = setdiff(ind_D,ind_B_fix);

KK_D = KK_nobc(ind_D_fix,ind_D_fix);
KK_bc= KK_nobc(ind_D_fix,ind_B_fix);

% % tic

[LK,UK,PK,QK] = lu(KK_D);

rr_Gauss = P_Gauss(:,1);

tolleranza_fix = SETTINGS.TOLL_FIXB;
N_ITERMAX_FIXB  = SETTINGS.N_ITERMAX_FIXB;

residuo_fix = 100;
N_ITER_FIX  = 1;

Rplot=[.9*min(meshData_loc.n(:,1)) 1.1*max(meshData_loc.n(:,1))];
Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
rgrid=linspace(Rplot(1),Rplot(2),200);
zgrid=linspace(Zplot(1),Zplot(2),200);
[RR,ZZ] = meshgrid(rgrid,zgrid);
% % % % 
% % figure; surf(RR,ZZ,griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk.Psi,RR,ZZ))
% % hold on
% % plot3(meshData_loc.n(ind_n_ONpla,1),meshData_loc.n(ind_n_ONpla,2),solk.Psi(ind_n_ONpla),'o')
% % plot3(meshData_loc.n(ind_n_ONpla,1),meshData_loc.n(ind_n_ONpla,2),Psi_fix_old(ind_n_ONpla),'o')


% %     figure
% % % %     edges=pdemesh(meshData_loc.n',meshData_loc.e_interface',[]);
% % % %     set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
% %     axis([Rplot Zplot]),  hold on,
% %         PSI_1 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk.Psi,RR,ZZ);
% %         contour(RR,ZZ,PSI_1,100);
% % % %     contour(RR,ZZ,PSI_2,100); colormap cool;



%%

while residuo_fix > tolleranza_fix && N_ITER_FIX < N_ITERMAX_FIXB
         
    
% %     [Source,solk_fix] = fun_calcSourceFEM_adaptive(meshData_loc,solk_fix,plaparameter,SETTINGS);
% %     lambda_fix = Ipla/sum(Source);
% %     Source = lambda_fix*Source;
        
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2(tri,solk_fix.Psi,N_order,P_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        [Psi_Gauss] = fun_evaluateFluxGaussPoints_v2_mex(tri,solk_fix.Psi,N_order,P_Gauss,n_Gauss,shape_functions);
    end
        
    Psi_k=Psi_Gauss;
    Psi_axis=solk_fix.Psi_axis;
    Psi_Bpla= solk.Psi_B;
    psi_bar=(Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux
    solk_fix.psi_bar=psi_bar;
    
    % for current density: only Gauss nodes inside plasma
    if SETTINGS.J_PARAMETRIZATION_TYPE == 1
        psi_bar(psi_bar>1) = NaN;
        psi_bar(psi_bar<0) = NaN;
        
        fdfn=interp1(psibar,FdF,psi_bar,'linear');  
        dpn=interp1(psibar,dP,psi_bar,'linear');   
        gg=(rr_Gauss.*dpn+fdfn./(mu0*rr_Gauss)); %computation of nodes plasma current density
        gg=gg/2/pi;
        
        gg(isnan(gg)) = 0;
        
        ind_t_selected = meshData_loc.ind_t_INpla;        

        if SETTINGS.RUN_MEX_ROUTINE == false
            [Source] = fun_calcSourceFEM_v2(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
            [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
        end

        Source(isnan(Source)) = 0;
        lambda_fix = Ipla/sum(Source);
        Source = lambda_fix*Source;
        
    elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
        tmp_gg = ((1-psi_bar.^alpha_M).^alpha_N);
        tmp_gg(imag(tmp_gg) ~= 0 ) = abs(tmp_gg(imag(tmp_gg) ~= 0)); % occhio!!!
        aa_r = (rr_Gauss*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss);
        gg = tmp_gg.*aa_r;
        
        ind_t_selected = meshData_loc.ind_t_INpla;
        
        if SETTINGS.RUN_MEX_ROUTINE == false
            [Source] = fun_calcSourceFEM_v2(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
            [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
        end
        
        lambda_fix = Ipla/sum(Source);
        Source = lambda_fix*Source;
        
    end
    
% %     F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),gg);
% % % % % %     F = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),psi_bar);
% % % %     % %         J_interp = zeros(meshData_loc.nn,1);
% %     J_interp = F(meshData_loc.n(:,1),meshData_loc.n(:,2));
% % % %     
% %     Jpla.faces=meshData_loc.t(:,1:3);
% %     Jpla.vertices=meshData_loc.n;
% %     Jpla.facevertexcdata=J_interp;
% %     
% %     figure
% %     axis([Rplot Zplot]),  hold on,
% %     hh=patch(Jpla,'facecolor','interp','edgecolor','none');
% %     axis equal;  hold on; colorbar vert; % title('j_\phi');
% %     axis([Rplot Zplot])
    
    % %         figure;  plot3(meshData_loc.nodes_Gauss_pla(qq,1),meshData_loc.nodes_Gauss_pla(qq,2),gg_tot(qq),'o')
    % %         figure;  plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),Source,'o')
    % %         figure;  plot3(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),gg_tot,'o')
    
    
    %
          
    Psi_D_fix = QK*(UK\(LK\(PK*(Source(ind_D_fix)*mu0 - KK_bc*Psi_B_fix))));
    
    % %     PSI_2 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk_fix.Psi,RR,ZZ);
% %     figure; surf(RR,ZZ,griddata(meshData_loc.n(ind_D_fix,1),meshData_loc.n(ind_D_fix,2),Psi_D_fix,RR,ZZ));

  
    Psi_fix = zeros(meshData_loc.nn,1);
    Psi_fix(ind_D_fix) = Psi_D_fix;
    Psi_fix(ind_B_fix) = Psi_B_fix;
        
    solk1_fix.Psi = Psi_fix;
    solk1_fix.lambda = lambda_fix;
    solk1_fix.lambda_star = mu0*lambda_fix;
    
% %     solk1_fix.Psi_B = solk_fix.Psi_B;
% %     solk1_fix.Separatrix = solk_fix.Separatrix;
% %     solk1_fix.Psi_B = solk_fix.Psi_B;
% %     solk1_fix.ind_t_Xp = solk_fix.ind_t_Xp;
    
    
    [solk1_fix] = fun_FluxGauss(meshData_loc,solk1_fix,SETTINGS);
    [solk1_fix] = fun_GradientGauss(meshData_loc,solk1_fix,SETTINGS);

    % update magnetic axis
    [solk1_fix]=fun_MagneticAxis(meshData_loc,solk1_fix,SETTINGS);

    solk1_fix.Psi_axis;

    residuo_fix = norm(Psi_fix - Psi_fix_old);
    
    solk_fix = solk1_fix;
    solk_fix.N_ITER_FIX = N_ITER_FIX;
    Psi_fix_old = Psi_fix;
    N_ITER_FIX = N_ITER_FIX+1;
    
% %     PSI_2 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk1_fix.Psi,RR,ZZ);
% % % %     figure; surf(RR,ZZ,PSI_2);
% %     
% %     figure
% %     Rplot=[.9*min(meshData_loc.n(:,1)) 1.1*max(meshData_loc.n(:,1))];
% %     Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
% %     rgrid=linspace(Rplot(1),Rplot(2),200);
% %     zgrid=linspace(Zplot(1),Zplot(2),200);
% %     contour(RR,ZZ,PSI_2,[solk.Psi_Bpla solk.Psi_Bpla],'r','LineWidth',2);
% %     [RR,ZZ] = meshgrid(rgrid,zgrid);
% %     axis([Rplot Zplot]),  hold on,
% %     contour(RR,ZZ,PSI_2,50);
% %     colorbar vert; colormap jet
% %     plot(meshData_loc.n(solk1_fix.ind_n_axis,1),meshData_loc.n(solk1_fix.ind_n_axis,2),'o')
% %     quiver(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2), ...
% %         solk1_fix.Grad_Gauss(:,1),solk1_fix.Grad_Gauss(:,2))
% %     axis equal
% % % %     
% %       pause
    
end
% % toc
% % aa = 0;
% %     Rplot=[.9*min(meshData_loc.n(:,1)) 1.1*max(meshData_loc.n(:,1))];
% %     Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
% %     rgrid=linspace(Rplot(1),Rplot(2),200); 
% %     zgrid=linspace(Zplot(1),Zplot(2),200);
% %     [RR,ZZ] = meshgrid(rgrid,zgrid);
% % 
% % PSI_2 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),solk_fix.Psi,RR,ZZ);
% % figure
% % hold on
% %  contour(RR,ZZ,PSI_2,50);
% %     colorbar vert; colormap jet
% %  contour(RR,ZZ,PSI_2,solk.Psi_B*[1 1]);
% % 
% % plot3(meshData_loc.n(ind_B_fix,1),meshData_loc.n(ind_B_fix,2),Psi_B_fix,'ro')
% % contour3(RR,ZZ,PSI_2,[solk.Psi_Bpla solk.Psi_Bpla],'r','LineWidth',2);
% % figure
% % surf(RR,ZZ,PSI_2);


% % Rplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,1))...
% %     max(meshData_loc.n(meshData_loc.ind_n_bc,1))]+[-0.3 0.3];
% % Zplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,2))...
% %     max(meshData_loc.n(meshData_loc.ind_n_bc,2))]+[-0.3 0.3];
% % 
% % rr=meshData_loc.n(:,1);
% % zz=meshData_loc.n(:,2);
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % 
% % gg_pla_int = scatteredInterpolant(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),gg);
% % gg_pla_nodes = gg_pla_int(meshData_loc.n(:,1),meshData_loc.n(:,2));
% % 
% % Jpla.faces=meshData_loc.t(:,[1:3 1]);
% % Jpla.vertices=meshData_loc.n;
% % Jpla.facevertexcdata=Source;
% % 
% % 
% % figure
% % % % edges=pdemesh(meshData.n',meshData.e',[]);
% % % % set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
% % hh=patch(Jpla,'facecolor','interp','edgecolor','none');
% % axis([Rplot Zplot]),  hold on,
% % axis equal;  hold on; colorbar vert; % title('j_\phi');
% % plot(solk_fix.Separatrix(:,1),solk_fix.Separatrix(:,2),'-r');
% % plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'-k')
% % colormap jet


% % PSI_2 = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),Psi_fix,RR,ZZ);
% % % % 
% % figure
% % % % edges=pdemesh(meshData.n',meshData.e',[]);
% % % % set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
% % Rplot=[.9*min(meshData_loc.n(:,1)) 1.1*max(meshData_loc.n(:,1))];
% % Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % axis([Rplot Zplot]),  hold on,
% % contour(RR,ZZ,PSI_2,100);
% % contour(RR,ZZ,PSI_2,[solk.Psi_Bpla solk.Psi_Bpla],'r','LineWidth',2);
% % colormap('cool'); colorbar vert;

% % figure
% % edges=pdemesh(meshData.n',meshData.e',[]);
% % set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
% % axis([Rplot Zplot]),  hold on,
% % JP.faces=meshData_loc.t(:,1:3);
% % JP.vertices=meshData_loc.n;
% % JP.facevertexcdata=Source./meshData_loc.Area_n;
% % hh=patch(JP,'facecolor','interp','edgecolor','none');
% % axis equal;  hold on; colorbar vert;  title('j_\phi'); colormap cool
% % contour(RR,ZZ,PSI_2,[solk.Psi_Bpla solk.Psi_Bpla],'r','LineWidth',2);
% % axis([Rplot Zplot])
