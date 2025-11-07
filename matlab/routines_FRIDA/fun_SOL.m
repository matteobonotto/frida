function [Separatrix]=fun_SOL(meshData,solk,SETTINGS)
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

meshData_SOL.t=[meshData.t(ind_global_SOL,1:3)';ones(1,size(meshData.t(ind_global_SOL,:),1))];
% %     triplot(meshData_SOL.t(1:3,:)',meshData.n(:,1),meshData.n(:,2))
Psi_interp=pdeInterpolant(meshData.n',meshData_SOL.t,solk.Psi);

% % vessel=meshData.n(meshData.vess,:);
FW_rz = meshData.n(meshData.ind_n_FW,:);
FW_rz_extended = repmat(FW_rz,3,1);

[FW_rz_sort,ind_FW_sort] = fun_ordinapunti(FW_rz);
FW_rz_sort_extended = repmat(FW_rz_sort,3,1);

dim_tri=((FW_rz(1:end-1,1)-FW_rz(2:end,1)).^2+...
    (FW_rz(1:end-1,2)-FW_rz(2:end,2)).^2).^0.5;
dim_tri=sum(dim_tri)/numel(dim_tri); %dimensione media triangoli

Rplot=[min(meshData.n(:,1))...
    max(meshData.n(:,1))]+[-0.3 0.3];
Zplot=[min(meshData.n(:,2))...
    max(meshData.n(:,2))]+[-0.3 0.3];
rgrid=linspace(Rplot(1),Rplot(2),200);
zgrid=linspace(Zplot(1),Zplot(2),200);


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
% % triplot(meshData_SOL.t(1:3,:)',meshData.n(:,1),meshData.n(:,2),'g')
% % contour(RR,ZZ,PSI,[min(solk.Psi_B) max(solk.Psi_B)],'k');
% % contour(RR,ZZ,PSI,[SOL(1) SOL(1)],'k');
% % % % contour(RR,ZZ,PSI,100,'k');
% % contour(RR,ZZ,PSI,[SOL(2) SOL(2)],'k');

%%

%%
if numel(solk.Psi_B)==1
    
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

% %         plot(retta(:,1),retta(:,2),'.k'); plot(retta(1,1),retta(1,2),'o')
        
% %                         plot(retta(:,1),retta(:,2),'.k'); plot(retta(1,1),retta(1,2),'o')
        uq_Cell = Psi_interp.evaluate(retta(:,1),retta(:,2));
        % %         ind_eval=find(isnan(uq_Cell)==0);
        ind_eval=find(~isnan(uq_Cell));
        uq_Cell=uq_Cell(ind_eval);
        retta=retta(ind_eval,:);
        Psi_star=solk.Psi_B;
        hh=1;
        xx=1:length(uq_Cell);
        xx_q=1:0.1:length(uq_Cell);
        uq_Cell=interp1(xx,uq_Cell,xx_q)';
        retta=interp1(xx,retta,xx_q);
        derivative=gradient(abs(uq_Cell-Psi_star));
        % %         figure
        % %         plot(abs(uq_Cell-Psi_star))
        % %         hold on
        % %         plot(derivative)
        
        if  max(sign(derivative))==min(sign(derivative))
            % %             while  max(sign(derivative))==min(sign(derivative))==1
            [~,Raggio_new]=cart2pol(retta(1,1)-solk.Axis_RR,...
                retta(1,2)-solk.Axis_ZZ);
            rr=linspace(Raggio_new,Raggio,100)';
            retta=[rr*cos(theta(jj)) rr*sin(theta(jj))];
            retta=[retta(:,1)+solk.Axis_RR ...
                retta(:,2)+solk.Axis_ZZ];
            % %                     plot(retta(:,1),retta(:,2),'.k')
            uq_Cell = Psi_interp.evaluate(retta(:,1),retta(:,2));
            % %             ind_eval=find(isnan(uq_Cell)==0);
            ind_eval=find(~isnan(uq_Cell));
            uq_Cell=uq_Cell(ind_eval);
            retta=retta(ind_eval,:);
            Psi_star=solk.Psi_B;
            hh=1;
            derivative=gradient(abs(Psi_star-uq_Cell));
            % %             end
        end
        
        while derivative(hh)<0 && hh<numel(derivative)
            hh=hh+1;
        end
        
        % %         abs(Psi_star-uq_Cell(hh))
        diff_psi=abs(Psi_star-uq_Cell);
        % %         jj
        
        if hh<length(diff_psi)
            if diff_psi(hh-1)<diff_psi(hh)
                hh=hh-1;
            elseif diff_psi(hh+1)<diff_psi(hh)
                hh=hh+1;
            end
        end
        
        % %         abs(Psi_star-uq_Cell(hh))
        % %         if abs(Psi_star-uq_Cell(hh))<0.01
        % %             ind=hh;
        % %         end
        ind=hh;
        Punti(jj,:)=retta(ind,:);
% %         plot(retta(:,1),retta(:,2),'.k')
% %         plot(Punti(jj,1),Punti(jj,2),'*r')
% %                 pause
% %         pause
        
    end
    
    xP=1:numel(Punti(:,1));
    xQ=linspace(1,numel(Punti(:,1)),n_points_Separatrix);
    Separatrix=spline(xP,Punti',xQ)';
% %     Separatrix=interp1(xP',Punti,xQ');
    % %     plot(Separatrix(:,1),Separatrix(:,2),'.-k')
    
elseif numel(solk.Psi_B)>1
    % PROXIMITY PARAMETER
    solk.plasma_radius=.5;
    dist_Xpoints=sqrt((solk.XP_RR(1)-solk.XP_RR(2))^2+...
        (solk.XP_ZZ(1)-solk.XP_ZZ(2))^2);
    solk.sigma=dist_Xpoints/solk.plasma_radius;
    
    if solk.sigma<.5
        [~,ind]=min(abs(solk.Psi_axis-solk.Psi_B));
        solk.Psi_B=solk.Psi_B(ind);
        solk.XP_RR=solk.XP_RR(ind);
        solk.XP_ZZ=solk.XP_ZZ(ind);
    end
    
    
    %% Multiple Xpoint
    
    Xpoints=[solk.XP_RR solk.XP_ZZ];
    
    for ii=1:numel(solk.Psi_B)
        [~,R]=cart2pol(Xpoints(ii,1)-solk.Axis_RR,...
            Xpoints(ii,2)-solk.Axis_ZZ);
        Raggio(ii,1)=R;
        TH(ii)=fun_myatan(Xpoints(ii,2)-solk.Axis_ZZ,...
            Xpoints(ii,1)-solk.Axis_RR);
    end
    
    [TH,ind]=sort(TH);
    Xpoints=Xpoints(ind,:);
    solk.Psi_B=solk.Psi_B(ind);
    
    % Iterative procedure to find separatrix in a piece-wise way
    Separatrix=[];
    
    for ii=1:numel(solk.Psi_B)
        
        if ii<numel(solk.Psi_B)
            TH_1=TH(ii);
            TH_2=TH(ii+1);
            XP_Point_1=Xpoints(ii,:);
            XP_Psi_1=solk.Psi_B(ii);
            XP_Point_2=Xpoints(ii+1,:);
            XP_Psi_2=solk.Psi_B(ii+1);
            n_angle=ceil(abs(TH_1-TH_2)/(pi/6));
            theta=linspace(TH_1,TH_2,n_angle);
        elseif ii==numel(solk.Psi_B)
            TH_1=TH(end)-2*pi;
            TH_2=TH(1);
            XP_Point_1=Xpoints(end,:);
            XP_Psi_1=solk.Psi_B(end);
            XP_Point_2=Xpoints(1,:);
            XP_Psi_2=solk.Psi_B(1);
            n_angle=ceil(abs(TH_1-TH_2)/(pi/6));
            theta=linspace(TH_1,TH_2,n_angle);
            % %             plot(retta(:,1),retta(:,2),'.k')
            % %             plot(Punti(jj,1),Punti(jj,2),'*r')
            % %             pause
        end
        
        nPunti=numel(theta);
        Delta_Psi=(XP_Psi_2-XP_Psi_1)/(nPunti);
        Punti=zeros(nPunti,2);
        Punti(1,:)=XP_Point_1;
        Punti(end,:)=XP_Point_2;
        
        for jj=2:nPunti-1
            rr=linspace(0,max(Raggio),500)';
            retta=[rr*cos(theta(jj)) rr*sin(theta(jj))];
            retta=[retta(:,1)+solk.Axis_RR ...
                retta(:,2)+solk.Axis_ZZ];
            % %                 plot(retta(:,1),retta(:,2),'.k')
            uq_Cell = Psi_interp.evaluate(retta(:,1),retta(:,2));
            ind_eval=find(isnan(uq_Cell)==0);
            uq_Cell=uq_Cell(ind_eval);
            retta=retta(ind_eval,:);
            Psi_star=XP_Psi_1+(jj-1)*Delta_Psi;
            hh=1; xx=1:length(uq_Cell);
            xx_q=1:0.2:length(uq_Cell);
            uq_Cell=interp1(xx,uq_Cell,xx_q)';
            retta=interp1(xx,retta,xx_q);
            derivative=gradient(abs(uq_Cell-Psi_star));
            if  max(sign(derivative))==min(sign(derivative))
                % %                 while  max(sign(derivative))==min(sign(derivative))==1
                [~,Raggio_new]=cart2pol(retta(1,1)-solk.Axis_RR,...
                    retta(1,2)-solk.Axis_ZZ);
                rr=linspace(Raggio_new,Raggio,100)';
                retta=[rr*cos(theta(jj)) rr*sin(theta(jj))];
                retta=[retta(:,1)+solk.Axis_RR ...
                    retta(:,2)+solk.Axis_ZZ];
                % %         plot(retta(:,1),retta(:,2),'.k')
                uq_Cell = Psi_interp.evaluate(retta(:,1),retta(:,2));
                ind_eval=find(isnan(uq_Cell)==0);
                uq_Cell=uq_Cell(ind_eval);
                retta=retta(ind_eval,:);
                Psi_star=solk.Psi_B;
                hh=1;
                derivative=gradient(abs(Psi_star-uq_Cell));
                % %                 end
            end
            
            while derivative(hh)<0 && hh<numel(derivative)
                hh=hh+1;
            end
            
            % %         abs(Psi_star-uq_Cell(hh))
            diff_psi=abs(Psi_star-uq_Cell);
            if diff_psi(hh-1)<diff_psi(hh)
                hh=hh-1;
            elseif diff_psi(hh+1)<diff_psi(hh)
                hh=hh+1;
            end
            
            % %         abs(Psi_star-uq_Cell(hh))
            % %         if abs(Psi_star-uq_Cell(hh))<0.01
            % %             ind=hh;
            % %         end
            ind=hh;
            Punti(jj,:)=retta(ind,:);
            % %             plot(retta(:,1),retta(:,2),'.k')
            % %             plot(Punti(jj,1),Punti(jj,2),'*r')
            % %             pause
            
        end
        xP=1:numel(Punti(:,1));
        xQ=linspace(1,numel(Punti(:,1)),n_points_Separatrix);
        Separatrix_jj=spline(xP,Punti',xQ)';
        % %         plot(Separatrix_jj(:,1),Separatrix_jj(:,2),'k.')
        Separatrix=[Separatrix; Separatrix_jj(2:end,:)];
        
    end
    
    Separatrix=[Separatrix; Separatrix(1,:)];
    
end



end
