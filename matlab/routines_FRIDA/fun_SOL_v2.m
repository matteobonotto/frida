function [Separatrix]=fun_SOL_v2(meshData,solk,SETTINGS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find the boundary ([r,z] curve) for a given plasma configuration.
%
%   If the plasma is Limiter/Xpoint with a single null, the routine returns
%   the iso-flux curve for the passing for the constrain point (Limiter
%   point or Xpoint).
%
%   If there are multiple Xpoints the procedure finds the curve that
%   connect these Xpoints, even if the flux value at these point are
%   different (the curve is no longher iso-flux curve).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('DATA.mat')
ind_t_Invess = meshData.ind_t_InFW;
ind_t_Invess = 1:meshData.nt;
index_n_vess = meshData.ind_n_InFW;

n_points_Separatrix = SETTINGS.N_POINTS_SEPARATRIX;

%% Define region where search for boundary

if ~any(strcmp('fac_tresh_SOL',fieldnames(SETTINGS)))
    fac_tresh_SOL = 3;
else
    fac_tresh_SOL = SETTINGS.fac_tresh_SOL;
end
SOL=[min(solk.Psi(index_n_vess)) max(solk.Psi_B)+(solk.Psi_axis-min(solk.Psi_B))/fac_tresh_SOL];
Psi_trianges=solk.Psi(meshData.t(ind_t_Invess,:));
Psi_centro_t=sum(Psi_trianges,2)/size(Psi_trianges,2);
ind_t_SOL=find(Psi_centro_t<max(SOL) & Psi_centro_t>min(SOL));

ind_global_SOL = ind_t_Invess(ind_t_SOL);%ind_t_SOL;

meshData_SOL.t = meshData.t(ind_global_SOL,1:3);


% % meshData_SOL.t = meshData.t;

% %     triplot(meshData_SOL.t(1:3,:)',meshData.n(:,1),meshData.n(:,2))

% % Psi_interp=pdeInterpolant(meshData.n',meshData_SOL.t(:,[1:6 1]).',solk.Psi);
Psi_interp=pdeInterpolant(meshData.n',meshData_SOL.t(:,[1:3 1]).',solk.Psi);



% % vessel=meshData.n(meshData.vess,:);
FW_rz = meshData.n(meshData.ind_B,:);
FW_rz_extended = repmat(FW_rz,3,1);

[FW_rz_sort,ind_FW_sort] = fun_ordinapunti(FW_rz);
FW_rz_sort_extended = repmat(FW_rz_sort,3,1);

dim_tri=((FW_rz(1:end-1,1)-FW_rz(2:end,1)).^2+...
    (FW_rz(1:end-1,2)-FW_rz(2:end,2)).^2).^0.5;
dim_tri=sum(dim_tri)/numel(dim_tri); %dimensione media triangoli

% % Rplot=[min(meshData.n(:,1))...
% %     max(meshData.n(:,1))]+[-0.3 0.3];
% % Zplot=[min(meshData.n(:,2))...
% %     max(meshData.n(:,2))]+[-0.3 0.3];
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% %
% %
% % figure
% % hold on; axis equal; colormap jet
% % hold on, xlabel('r [m]'), ylabel('z [m]');
% % % % edges=pdemesh(meshData.n',meshData.e',[]);
% % % % set(edges,'color','k')
% % % % set(edges,'linewidth',1)
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % rr=meshData.n(:,1);
% % zz=meshData.n(:,2);
% % PSI = griddata(rr,zz,solk.Psi,RR,ZZ);
% % contour(RR,ZZ,PSI,100);
% % triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2),'r')
% % triplot(meshData_SOL.t,meshData.n(:,1),meshData.n(:,2),'g')
% % contour(RR,ZZ,PSI,[min(solk.Psi_B) max(solk.Psi_B)],'k', 'LineWidth', 2);

%%

%% Limiter or single Xpoint

Xpoints=[solk.XP_RR solk.XP_ZZ];
[~,Raggio_ves]=cart2pol(FW_rz(:,1)-solk.Axis_RR,FW_rz(:,2)-solk.Axis_ZZ);
TH_ves=zeros(size(FW_rz,1),1);

for hh=1:size(FW_rz,1)
    th_ves=fun_myatan(FW_rz(hh,2)-solk.Axis_ZZ,FW_rz(hh,1)-solk.Axis_RR);
    TH_ves(hh,1)=th_ves;
end

TH_ves_extended = [TH_ves-2*pi; TH_ves; TH_ves+2*pi];
[TH_ves_extended,ind_TH_ves_extended_sort] = sort(TH_ves_extended);
FW_rz_sort_extended = FW_rz_extended(ind_TH_ves_extended_sort,:);

TH=fun_myatan(Xpoints(2)-solk.Axis_ZZ,Xpoints(1)-solk.Axis_RR);
theta=(TH(1):pi/20:TH(1)+2*pi);
nPunti=numel(theta);
Punti=zeros(nPunti,2);
Punti(1,:)=Xpoints;
Punti(end,:)=Xpoints;


for jj=2:nPunti-1
    if theta(jj)>2*pi
        theta(jj)=theta(jj)-2*pi;
    end
    
    FW_star = interp1(TH_ves_extended,FW_rz_sort_extended,theta(jj));
    
    if any(isnan(FW_star)) || any(isempty(FW_star))
        [~,ind]=min(abs(TH_ves-theta(jj)));
        Raggio=Raggio_ves(ind);
    else
        Raggio = norm(FW_star - [solk.Axis_RR solk.Axis_ZZ]);
    end
    
    rr=linspace(.01*Raggio,Raggio,150)';
    retta=[rr*cos(theta(jj)) rr*sin(theta(jj))];
    retta=[retta(:,1)+solk.Axis_RR,retta(:,2)+solk.Axis_ZZ];
    
    Psi_star=solk.Psi_B;
    
    uq_Cell = Psi_interp.evaluate(retta(:,1),retta(:,2));
    s_coo = [0 ;cumsum(fun_vecnorm(retta(2:end,:)-retta(1:end-1,:)))];
    ind_sel = ~isnan(uq_Cell);
    
    % %         figure
    % %         plot(s_coo,abs(uq_Cell-Psi_star),'o');
    % %         hold on
    % %         plot(s_coo(ind_sel),abs(uq_Cell(ind_sel)-Psi_star),'*');
    
    
    [s_coo_star,~] = fun_MacFindX(s_coo(ind_sel),uq_Cell(ind_sel),Psi_star);
    % %     if length(s_coo_star) > 1; s_coo_star= max(s_coo_star); end
    
    if isempty(s_coo_star)
        rr=linspace(.01*Raggio,Raggio,2*150)';
        retta=[rr*cos(theta(jj)) rr*sin(theta(jj))];
        retta=[retta(:,1)+solk.Axis_RR,retta(:,2)+solk.Axis_ZZ];
        
        Psi_star=solk.Psi_B;
        uq_Cell = Psi_interp.evaluate(retta(:,1),retta(:,2));
        s_coo = [0 ;cumsum(fun_vecnorm(retta(2:end,:)-retta(1:end-1,:)))];
        
        [s_coo_star,~] = fun_MacFindX(s_coo(ind_sel),uq_Cell(ind_sel),Psi_star);
        % %         if length(s_coo_star) > 1; s_coo_star= max(s_coo_star); end
        
    end
    
    if numel(s_coo_star) > 1
        s_coo_star = min(s_coo_star);
    end
    
    % %     jj
    try
        Punti(jj,:)=[interp1(s_coo,retta(:,1),s_coo_star) interp1(s_coo,retta(:,2),s_coo_star)];
    catch
        aa = 0;
    end
    % %     plot(retta(:,1),retta(:,2),'.k'); plot(Punti(jj,1),Punti(jj,2),'*r')
    % %                 pause
    % %         pause
    
end

xP=1:numel(Punti(:,1));
xQ=linspace(1,numel(Punti(:,1)),n_points_Separatrix);
Separatrix=spline(xP,Punti',xQ)';
% %     Separatrix=interp1(xP',Punti,xQ');
% %     plot(Separatrix(:,1),Separatrix(:,2),'.-k')






end









