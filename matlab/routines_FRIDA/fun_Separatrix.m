function [solk]=fun_Separatrix(meshData,solk,SETTINGS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find Plasma configuration (Limiter/Xpoint), multiple Xpoints are
%   supported.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% % ind_t_Invess=index.ind_t_Invess;
% % pp=meshData.n';

%% Nodes inside the vessel and on vessel


Mat_ww_Gauss = reshape(meshData.ww_Gauss_pla,meshData.n_Gauss,meshData.nt)';

Area_tri = sum(Mat_ww_Gauss,2); 

Area_tri_media = sum(Area_tri)/meshData.nt;

dim_tri = sqrt(2*Area_tri_media);

vessel=meshData.n(meshData.ind_n_FW,:);
% % dim_tri=((vessel(1:end-1,1)-vessel(2:end,1)).^2+...
% %     (vessel(1:end-1,2)-vessel(2:end,2)).^2).^0.5;
% % dim_tri=sum(dim_tri)/numel(dim_tri); %dimensione media triangoli

%% Remove magnetic axis from search area
if any(strcmp('ell_par',fieldnames(solk)))==1
    % If first iteration the equivalent plasma (elliptical) requires a
    % bigger area to be removed
    r0=solk.ell_par(1);
    z0=solk.ell_par(2);
    a=solk.ell_par(3);
    b=solk.ell_par(4);
    Axis=[r0 z0];
    disk_rem.RR=r0+1.1*a*cos(linspace(0,2*pi,100)');
    disk_rem.ZZ=z0+1.1*b*sin(linspace(0,2*pi,100)');
    % %     plot(disk_rem.RR,disk_rem.ZZ,'k.')
else
    
    Axis=[solk.Axis_RR solk.Axis_ZZ];
    % %     plot(Axis(1),Axis(2),'k*')
    disk_rem.RR=Axis(1)+5*dim_tri*cos(linspace(0,2*pi,100));
    disk_rem.ZZ=Axis(2)+5*dim_tri*sin(linspace(0,2*pi,100));
    % %     plot(disk_rem.RR,disk_rem.ZZ,'k.')
    
end

index_t = meshData.ind_t_InFW;
ind_remove = find(inpolygon(meshData.c_t(:,1),meshData.c_t(:,2),disk_rem.RR,disk_rem.ZZ));
ind_geo_ok=setdiff(index_t,ind_remove)';

% % triplot(meshData.t(ind_geo_ok,1:3),meshData.n(:,1),meshData.n(:,2),'g')

%%
Rplot=[min(meshData.n(meshData.ind_n_bc,1))...
    max(meshData.n(meshData.ind_n_bc,1))]+[-0.2 0.3];
Zplot=[min(meshData.n(meshData.ind_n_bc,2))...
    max(meshData.n(meshData.ind_n_bc,2))]+[-0.3 0.3];

rr=meshData.n(:,1);
zz=meshData.n(:,2);

rgrid=linspace(Rplot(1),Rplot(2),200);
zgrid=linspace(Zplot(1),Zplot(2),200);

% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % 
% %     figure
% %     hold on;    grid on;   axis equal; colormap jet; box on
% %     hold on, xlabel('r [m]'), ylabel('z [m]');
% %     axis([Rplot Zplot])
% %     contour(RR,ZZ,PSI,50);
    
% %     quiver(P_Gauss(:,1),P_Gauss(:,2),Grad_Gauss(:,1),Grad_Gauss(:,2),2)



%% Find X-Point(s)
% Iterative procedure to find all X-Point(s) inside plasma region


% % Psi_nodes        = solk.Psi;
% % tri              = meshData.t;
% % nodes            = meshData.n;
% % shape_functions  = meshData.shape_functions;
% % N_order          = meshData.shape_functions_N_order;
% % n_Gauss          = meshData.n_Gauss;
% % P_Gauss          = meshData.nodes_Gauss_pla;
% % 
% % if SETTINGS.RUN_MEX_ROUTINE == false
% %     [Grad_Gauss] = fun_calcGradientGauss_FEM(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
% % elseif SETTINGS.RUN_MEX_ROUTINE == true
% %     [Grad_Gauss] = fun_calcGradientGauss_FEM_mex(tri,Psi_nodes,P_Gauss,shape_functions,N_order,n_Gauss);
% % end


N_order = meshData.shape_functions_N_order;
Psi_Gauss  = solk.Psi_Gauss;
Grad_Gauss = solk.Grad_Gauss;
P_Gauss    = meshData.nodes_Gauss_pla;
n_Gauss   = meshData.n_Gauss;

XP_Psi = [];
XP_RR  = [];
XP_ZZ  = [];
ind_t_Xp = [];

min_grad_tot = [];

distanza = 100;
disk_circ_remove = 5;

IND_Gauss_all = reshape(1:meshData.nt*n_Gauss,n_Gauss,meshData.nt)';

min_Grad=0;
iter=0;
while min_Grad < 0.1 %|| iter < 2 % min_Grad < 0.05
    % New Xpoint    
    IND_Gauss_ok = IND_Gauss_all(ind_geo_ok,:);
    ind_Gauss_ok = IND_Gauss_ok(:);
    
    [min_grad,ind_G_min_ok]=min(abs(Grad_Gauss(ind_Gauss_ok,1)+1i*Grad_Gauss(ind_Gauss_ok,2)));
    [max_grad,~]=max(abs(Grad_Gauss(ind_Gauss_ok,1)+1i*Grad_Gauss(ind_Gauss_ok,2)));
    
    ind_G_min = ind_Gauss_ok(ind_G_min_ok);
% %     plot(P_Gauss(ind_G_min,1), P_Gauss(ind_G_min,2), 'o')
    
    min_grad = min_grad/max_grad; % normalizing gradient
    
    ind_t_min = ceil(ind_G_min/n_Gauss);
    if ismember(ind_G_min,IND_Gauss_all(ind_t_min,:)) == 0; error('wrong selection of temporary Xp'); end
    
% %     triplot(meshData.t(ind_t_min,1:3),meshData.n(:,1),meshData.n(:,2),'r', 'linewidth', 2)

    XPinterp.Psi = Psi_Gauss(ind_G_min);
    XPinterp.RR = P_Gauss(ind_G_min,1);
    XPinterp.ZZ = P_Gauss(ind_G_min,2);
    min_Grad = min_grad;
        
    % Interpolation
    dist_vess=min(((meshData.c_t(ind_t_min,1)-vessel(:,1)).^2+(meshData.c_t(ind_t_min,2)-vessel(:,2)).^2).^0.5);
% %     if dist_vess>2*dim_tri && min_grad<0.2
    if  min_grad<0.2

        if min(((XPinterp.RR-XP_RR).^2+(XPinterp.ZZ-XP_ZZ).^2).^0.5)<4*dim_tri
            min_Grad=100;
        
        elseif ((Axis(1)-XPinterp.RR)^2+(Axis(2)-XPinterp.ZZ)^2)^0.5<4*dim_tri
            min_Grad=100;
            
        end
    else
        min_Grad=100;
    end
    
    min_grad_tot = [min_grad_tot; min_Grad];
    
    if norm(Axis - [XPinterp.RR XPinterp.ZZ]) < 10*dim_tri
        min_Grad = 100;
    end
    
    % Update Xpoint vector if it is true Xpoint
    if iter > 0
        distanza = sqrt((XP_RR(1)-XPinterp.RR)^2 + (XP_ZZ(1)-XPinterp.ZZ)^2);
    end
    
    % %     if min_Grad < 0.05 && min_Grad/min_grad_tot(1) < 20 && distanza>5*disk_circ_remove*dim_tri
    if min_Grad < 0.1 && distanza>5*disk_circ_remove*dim_tri
        
        XP_Psi=[XP_Psi; XPinterp.Psi];
        XP_RR=[XP_RR; XPinterp.RR];
        XP_ZZ=[XP_ZZ; XPinterp.ZZ];
        ind_t_Xp = [ind_t_Xp; ind_t_min];
        
        Pcirc=[XPinterp.RR XPinterp.ZZ];
        
        % Update ind_geo_ok tu find next Xpoint
        disk_rem.RR=Pcirc(1)+disk_circ_remove*dim_tri*cos(linspace(0,2*pi,100));
        disk_rem.ZZ=Pcirc(2)+disk_circ_remove*dim_tri*sin(linspace(0,2*pi,100));
        
        % %         plot(disk_rem.RR,disk_rem.ZZ,'k.')
        
        % %         ind_remove=(inpolygon(meshData.n(index_n,1),meshData.n(index_n,2),disk_rem.RR,disk_rem.ZZ));
        % %         ind_remove = index_n(ind_remove);
        % %         ind_geo_ok=setdiff(ind_geo_ok,ind_remove);
        
        ind_remove = (inpolygon(meshData.c_t(ind_geo_ok,1),meshData.c_t(ind_geo_ok,2),disk_rem.RR,disk_rem.ZZ));
        ind_remove = ind_geo_ok(ind_remove);
        ind_geo_ok=setdiff(ind_geo_ok,ind_remove)';
        
        % %         triplot(meshData.t(ind_remove,1:3),meshData.n(:,1),meshData.n(:,2),'g')
        % %         triplot(meshData.t(ind_geo_ok,1:3),meshData.n(:,1),meshData.n(:,2),'g')
        
        
    else
        min_Grad=100;
        
    end
    
    iter=iter+1;
    
    % %     plot(disk_rem.RR,disk_rem.ZZ,'r-')
    % %     triplot(meshData.t(ind_remove,1:3),meshData.n(:,1),meshData.n(:,2),'color',rand(3,1))
    % %     triplot(meshData.t(ind_geo_ok,1:3),meshData.n(:,1),meshData.n(:,2),'color',rand(3,1))
end

%% Limiter (if is not diverted)

% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % 
% %     figure
% %     hold on;    grid on;   axis equal; colormap jet; box on
% %     hold on, xlabel('r [m]'), ylabel('z [m]');
% %     axis([Rplot Zplot])
% %     contour(RR,ZZ,PSI,50);
    

% % curve = meshData.n(meshData.ind_n_FW,:);
% % psi_curve = solk.Psi(meshData.ind_n_FW);
% % grad_curve = solk.Grad_nodes(meshData.ind_n_FW,:);
% % 
% % if ~isfield(meshData, 'vers_n')
% %     
% %     [gamma_sort,ind_1] = fun_ordinapunti(meshData.n(meshData.ind_n_FW,:));
% %     FW_new = gamma_sort;
% %     
% %     FW_mid = .5*(FW_new + FW_new([2:end 1],:));
% %     
% %     vers_t = FW_mid - FW_mid([2:end 1],:);
% %     vers_n = [-vers_t(:,2) vers_t(:,1)];
% %     
% %     vers_n_r = scatteredInterpolant(FW_new(:,1),FW_new(:,2),vers_n(:,1));
% %     vers_n_z = scatteredInterpolant(FW_new(:,1),FW_new(:,2),vers_n(:,2));
% %     
% %     vers_n = [vers_n_r(meshData.n(meshData.ind_n_FW,:)) vers_n_z(meshData.n(meshData.ind_n_FW,:))];
% %     meshData.vers_n = vers_n./fun_vecnorm(vers_n);
% % 
% %     % %     figure; quiver(meshData.n(meshData.ind_n_FW,1),meshData.n(meshData.ind_n_FW,2),vers_n(:,1),vers_n(:,2))
% %     
% % end

if ~isfield(meshData, 'vers_n')
    
    [gamma_sort,~] = fun_ordinapunti(meshData.n(meshData.ind_n_FW,:));
    FW_new = gamma_sort;
    
    FW_mid = .5*(FW_new + FW_new([2:end 1],:));
    
    vers_t = FW_mid - FW_mid([2:end 1],:);
    vers_n = [-vers_t(:,2) vers_t(:,1)];
    
    vers_n_r = scatteredInterpolant(FW_new(:,1),FW_new(:,2),vers_n(:,1));
    vers_n_z = scatteredInterpolant(FW_new(:,1),FW_new(:,2),vers_n(:,2));
    
    vers_n = [vers_n_r(meshData.n(meshData.ind_n_FW,:)) vers_n_z(meshData.n(meshData.ind_n_FW,:))];
    meshData.vers_n = vers_n./fun_vecnorm(vers_n);

    % %     figure; quiver(meshData.n(meshData.ind_n_FW,1),meshData.n(meshData.ind_n_FW,2),vers_n(:,1),vers_n(:,2))
    
end

qq = unique(meshData.t(:,4:end));
ind_n_FW_1st = setdiff(meshData.ind_n_FW,qq);

ind_sel_n_FW_1st = fun_vec_find(ind_n_FW_1st,meshData.ind_n_FW);

curve = meshData.n(meshData.ind_n_FW(ind_sel_n_FW_1st),:);
psi_curve = solk.Psi(meshData.ind_n_FW(ind_sel_n_FW_1st));
grad_curve = solk.Grad_nodes(meshData.ind_n_FW(ind_sel_n_FW_1st),:);
vers_n_1st = meshData.vers_n(ind_sel_n_FW_1st,:);


% % aa = tic;
% % for ii = 1:1000
% % [Psi_Lim,RZ_Lim] = fun_find_limiter_point_v3(curve,meshData.vers_n,psi_curve,grad_curve);
% % end
% % bb = toc(aa)/ii


% % [Psi_Lim,RZ_Lim] = fun_find_limiter_point_v2(curve,meshData.vers_n,psi_curve,grad_curve);
% % Lim.RR = RZ_Lim(1);
% % Lim.ZZ = RZ_Lim(2);
% % Lim.Psi = Psi_Lim;

% % load data_temp_debug.mat

[Psi_Lim,RZ_Lim] = fun_find_limiter_point_v3(curve,vers_n_1st,psi_curve,grad_curve);
Lim.RR = RZ_Lim(1);
Lim.ZZ = RZ_Lim(2);
Lim.Psi = Psi_Lim;

dist = meshData.c_t(meshData.ind_t_InFW,:) - RZ_Lim;
dist = sqrt(dist(:,1).^2 + dist(:,2).^2);
[~,ind_t_Lim] = min(dist);
ind_t_Lim = meshData.ind_t_InFW(ind_t_Lim);



% % Rplot=[min(meshData.n(meshData.ind_n_bc,1))...
% %     max(meshData.n(meshData.ind_n_bc,1))]+[-0.2 0.3];
% % Zplot=[min(meshData.n(meshData.ind_n_bc,2))...
% %     max(meshData.n(meshData.ind_n_bc,2))]+[-0.3 0.3];
% % 
% % rr=meshData.n(:,1);
% % zz=meshData.n(:,2);
% % 
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % 
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % 
% %     figure
% %     hold on;    grid on;   axis equal; colormap jet; box on
% %     hold on, xlabel('r [m]'), ylabel('z [m]');
% %     axis([Rplot Zplot])
% %     contour(RR,ZZ,PSI,50);
% %     contour(RR,ZZ,PSI,[Psi_Lim Psi_Lim], 'g', 'LineWidth',2);
% %     plot(RZ_Lim(1),RZ_Lim(2),'o');
% %     
% % aa = 0;
% % 
% % triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2),'b')
% % triplot(meshData.t(ind_t_Lim,1:3),meshData.n(:,1),meshData.n(:,2),'c', 'LineWidth',2);
% % plot(Lim.RR,Lim.ZZ,'or','LineWidth',3)
% % contour(RR,ZZ,PSI,Lim.Psi*[1 1], 'r', 'LineWidth',2);
% % contour(RR,ZZ,PSI,Psi_Lim*[1 1], 'r', 'LineWidth',2);
% % triplot(meshData.t(ind_geo_ok,1:3),meshData.n(:,1),meshData.n(:,2),'g')
% % plot(curve(:,1),curve(:,2),'o')

%%
% % Rplot=[min(meshData.n(meshData.ind_n_FW,1))...
% %     max(meshData.n(meshData.ind_n_FW,1))]+[-0.3 0.3];
% % Zplot=[min(meshData.n(meshData.ind_n_FW,2))...
% %     max(meshData.n(meshData.ind_n_FW,2))]+[-0.3 0.3];
% % 
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % 
% % figure
% % hold on;    grid on;   axis equal; colormap jet
% % axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % rr=meshData.n(:,1);
% % zz=meshData.n(:,2);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % contour(RR,ZZ,PSI,100);
% % contour(RR,ZZ,PSI,[min(XP_Psi) max(XP_Psi)],'k');
% % contour(RR,ZZ,PSI,[Lim.Psi Lim.Psi],'r');
% % plot(meshData.n(meshData.ind_n_FW,1),meshData.n(meshData.ind_n_FW,2),'.k')

%% Find configuration

% Chose only one Xpoint: to be removed
% % if numel(XP_Psi) > 1
% %     
% %     [XP_Psi,ind_max] = max(XP_Psi);
% %     
% %     XP_RR = XP_RR(ind_max);
% %     XP_ZZ = XP_ZZ(ind_max);
% %     
% % end


% % Rplot=[min(meshData.n(meshData.ind_n_bc,1))...
% %     max(meshData.n(meshData.ind_n_bc,1))]+[-0.2 0.3];
% % Zplot=[min(meshData.n(meshData.ind_n_bc,2))...
% %     max(meshData.n(meshData.ind_n_bc,2))]+[-0.3 0.3];
% % 
% % rr=meshData.n(:,1);
% % zz=meshData.n(:,2);
% % 
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % 
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % 
% % figure
% % edges=pdemesh(meshData.n',meshData.e_interface',[]);
% % set(edges,'color','k')
% % set(edges,'linewidth',1)
% % hold on;    grid on;   axis equal; colormap jet; box on
% % axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
% % contour(RR,ZZ,PSI,50);
% % contour(RR,ZZ,PSI,[Lim.Psi Lim.Psi],'r', 'LineWidth',2);
% % contour(RR,ZZ,PSI,[XP_Psi XP_Psi],'b', 'LineWidth',2);


if Lim.Psi>=max(XP_Psi)
    
    solk.CONFIG='Limiter';
    
    solk.Psi_B=Lim.Psi;
    solk.XP_RR=Lim.RR;
    solk.XP_ZZ=Lim.ZZ;
    solk.ind_t_Xp=ind_t_Lim;
    
elseif isempty(XP_Psi)
    
    solk.CONFIG='Limiter';
    
    solk.Psi_B=Lim.Psi;
    solk.XP_RR=Lim.RR;
    solk.XP_ZZ=Lim.ZZ;
    solk.ind_t_Xp=ind_t_Lim;
    
else
    
    solk.CONFIG='Diverted';
    
    [~,ind] = max(XP_Psi);
    
    solk.Psi_B=XP_Psi(ind);
    solk.XP_RR=XP_RR(ind);
    solk.XP_ZZ=XP_ZZ(ind);
    solk.ind_t_Xp=ind_t_Xp(ind);

end

% % plot(solk.XP_RR,solk.XP_ZZ,'ro')

% Se 2 nulli sono ad un flusso molto diverso prendo quello pi� distante dal
% vessel
if abs(100*(min(solk.Psi_B)-max(solk.Psi_B))/max(solk.Psi_B))>1.5
    % %     for ii=1:numel(solk.Psi_B)
    % %         dist=((vessel(:,1)-solk.XP_RR(ii)).^2+(vessel(:,2)-solk.XP_ZZ(ii)).^2).^0.5;
    % %         distanza(ii,1)=min(dist);
    % %     end
    % %     [~,ind_true]=max(distanza);
    % %     solk.Psi_B=solk.Psi_B(ind_true);
    % %     solk.XP_RR=solk.XP_RR(ind_true);
    % %     solk.XP_ZZ=solk.XP_ZZ(ind_true);
    
    [Psi_sort,ind_sort]=sort(solk.Psi_B);
    solk.Psi_B=max(Psi_sort(end));
    solk.XP_RR=solk.XP_RR(ind_sort(end));
    solk.XP_ZZ=solk.XP_ZZ(ind_sort(end));
    
end

solk.Psi_Bpla=solk.Psi_B;

solk.ind_n_XP = zeros(length(solk.Psi_B),1);

for ii = 1:length(solk.Psi_B)
    
    dist_ii = sqrt((meshData.n(:,1)-solk.XP_RR(ii)).^2 + ...
        (meshData.n(:,2)-solk.XP_ZZ(ii)).^2);
    
    [~,bb] = min(dist_ii);
    
    % %     plot(meshData.n(bb,1),meshData.n(bb,2),'ro')
    solk.ind_n_XP(ii) = bb;
    
end


%% find weigths for XP/Limiter

t_inFW = meshData.t(meshData.ind_t_InFW,:);

if N_order == 2
    % %     switch solk.CONFIG
    % %         case 'Diverted'
            
            % %             triplot(meshData.t(solk.ind_t_Xp,1:3),meshData.n(:,1),meshData.n(:,2),'r', 'linewidth', 2)
            
            ind_tri = solk.ind_t_Xp;
            tri_XP = meshData.t(ind_tri,:);
            
            rr_XP = solk.XP_RR;
            zz_XP = solk.XP_ZZ;
            
            Psi_ii = solk.Psi(tri_XP);
            
            P1_ii = meshData.n(tri_XP(1),:);
            P2_ii = meshData.n(tri_XP(2),:);
            P3_ii = meshData.n(tri_XP(3),:);
            
            coeffs = meshData.shape_functions(ind_tri,:);
            coeffs = reshape(coeffs,3*N_order,3*N_order)';
            
            aa = coeffs(:,1);
            bb = coeffs(:,2);
            cc = coeffs(:,3);
            dd = coeffs(:,4);
            ee = coeffs(:,5);
            ff = coeffs(:,6);
            
            delta_b_weights = (aa*rr_XP^2 + bb*zz_XP^2 + ...
                cc*rr_XP*zz_XP + dd*rr_XP + ee*zz_XP + ff)';
            
            % %             abs(delta_b_weights*Psi_ii - solk.Psi_B)
            if abs(delta_b_weights*Psi_ii - solk.Psi_B)/solk.Psi_B > 1e-2
                warning('weights for Xpoint computed inaccurately');
                save data_debug
                error('check delta_b_weights')
            end
            
            solk.delta_b_weights = delta_b_weights;
            solk.ind_n_XP = tri_XP;
            
            
            % %         case 'Limiter'
            % %
            % % % %             plot(meshData.n(solk.ind_n_XP,1),meshData.n(solk.ind_n_XP,2),'o')
            % %
            % %             solk.delta_b_weights = 1;
            % %
            % %     end

    
end



%% Find Separatrix

% % Rplot=[min(meshData.n(meshData.ind_n_bc,1))...
% %     max(meshData.n(meshData.ind_n_bc,1))]+[-0.2 0.3];
% % Zplot=[min(meshData.n(meshData.ind_n_bc,2))...
% %     max(meshData.n(meshData.ind_n_bc,2))]+[-0.3 0.3];
% % 
% % rr=meshData.n(:,1);
% % zz=meshData.n(:,2);
% % 
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % 
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % 
% %     figure
% %     hold on;    grid on;   axis equal; colormap jet; box on
% %     hold on, xlabel('r [m]'), ylabel('z [m]');
% %     axis([Rplot Zplot])
% %     contour(RR,ZZ,PSI,50);
% %     contour(RR,ZZ,PSI,[solk.Psi_B solk.Psi_B], 'g', 'LineWidth',2);


% % [Separatrix]=fun_SOL(meshData,solk,SETTINGS);
[Separatrix]=fun_SOL_v2(meshData,solk,SETTINGS);
solk.Separatrix=Separatrix; % % plot(Separatrix(:,1),Separatrix(:,2),'.-k')


% % 
% % aa = tic;
% % for ii = 1:100
% % [Separatrix]=fun_SOL_v2(meshData,solk,SETTINGS);
% % end
% % bb = toc(aa)/100




