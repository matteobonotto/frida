function [solk0] = fun_firstfluxmap(solk0,meshData,meshData_loc,plaparameter,SETTINGS)


Rplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,1))...
    max(meshData_loc.n(meshData_loc.ind_n_bc,1))]+[-0.3 0.3];
Zplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,2))...
    max(meshData_loc.n(meshData_loc.ind_n_bc,2))]+[-0.3 0.3];



%%
if SETTINGS.INITIAL_PLASMA_MODEL == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % %     r0=plaparameter.Centroid(1)+.2;
    % %     z0=plaparameter.Centroid(2)-.02;
    r0=plaparameter.Centroid(1);
    z0=plaparameter.Centroid(2);
    clear source point
    % %     point.RR=nodes_flux(:,1);     point.ZZ=nodes_flux(:,2);
% %     point.RR=meshData_loc.n(:,1);     point.ZZ=meshData_loc.n(:,2);
% %     source.R=r0;
% %     source.Z=z0;
% %     source.current=ones(size(source.R));
% %     source.type='ring';
% %     [GG_p_n]=fun_Field_Loop_flux(source,point);
% %     GG_p_n.psi = GG_p_n.psi/2/pi;
% %     Psi_n_pla = GG_p_n.psi*(plaparameter.Ipla);
% %     Psi_n=solk0.Psi_MS+Psi_n_pla;% + 0.3*plaparameter.Ipla);%Iell;
    
    solk0.Riphi=r0;
    solk0.Ziphi=z0;
    solk0.Iphi=plaparameter.Ipla;

    clear source point
    Points =meshData_loc.n;
    source = [r0 z0] ;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        
        Psi_n_pla=fun_Green_Flux_Loop(source,Points,plaparameter.Ipla);
        Psi_n_pla = Psi_n_pla/2/pi;
        
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        
        R_points = Points(:,1);
        Z_points = Points(:,2);
        npt_point = length(R_points);
        
        OPT_PARALLEL = 1;
        n_threads = 24;
        
        R_source = source(:,1);
        Z_source = source(:,2);
        npt_source = length(R_source);
        I_source = double(plaparameter.Ipla);
        
        psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
        Psi_n_pla = psi/2/pi;
        
    end
    
    Psi_n = solk0.Psi_MS + Psi_n_pla;% + 0.3*plaparameter.Ipla);%Iell;
    
    
    % %     figure
    % %     edges=pdemesh(meshData.n',meshData.e_interface',[]);
    % %     set(edges,'color','k')
    % %     set(edges,'linewidth',1)
    % %     rr=meshData_loc.n(:,1);
    % %     zz=meshData_loc.n(:,2);
    % %     axis equal;
    % %     axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
    % %     rgrid=linspace(Rplot(1),Rplot(2),200);
    % %     zgrid=linspace(Zplot(1),Zplot(2),200);
    % %     [RR,ZZ] = meshgrid(rgrid,zgrid);
    % %     PSI = griddata(rr,zz,Psi_n_pla,RR,ZZ) - griddata(rr,-zz,Psi_n_pla,RR,ZZ);
    % %     contour(RR,ZZ,PSI,100)
    
    
elseif SETTINGS.INITIAL_PLASMA_MODEL == 2
    
    % Pla (Ne equivalent currents)
    Ne=40;
    r0=plaparameter.Centroid(1);
    z0=plaparameter.Centroid(2);
    % %     a=0.15*(max(meshData.n(meshData.ind_n_FW,1))-min(meshData.n(meshData.ind_n_FW,1)));
    % %     b=0.15*(max(meshData.n(meshData.ind_n_FW,2))-min(meshData.n(meshData.ind_n_FW,2)));
    a=0.15*(max(meshData_loc.n(meshData_loc.ind_B,1))-min(meshData_loc.n(meshData_loc.ind_B,1)));
    b=0.15*(max(meshData_loc.n(meshData_loc.ind_B,2))-min(meshData_loc.n(meshData_loc.ind_B,2)));
    theta=linspace(0,2*pi,Ne+1)';
    theta = theta(1:end-1);
    ell.RR=r0+a*cos(theta);
    ell.ZZ=z0+b*sin(theta);
    
    dist = min(sqrt((r0 - meshData_loc.n(meshData_loc.ind_n_FW,1)).^2 ...
        + (z0 - meshData_loc.n(meshData_loc.ind_n_FW,2)).^2));
    
    if .5*min([a b]) > dist
        a = .2*dist;
        b = .2*dist;
        
        ell.RR=r0+a*cos(theta);
        ell.ZZ=z0+b*sin(theta);
        
    end
    
    solk0.ell_par=[r0 z0 a b];
    solk0.Riphi=ell.RR;
    solk0.Ziphi=ell.ZZ;
    
    
    % Match moments of current density
    I0=plaparameter.Ipla/Ne;
    Is=2*plaparameter.Ipla*(plaparameter.Centroid(2)-z0)/(b*Ne);
    Ic=plaparameter.Ipla*(plaparameter.Centroid(1)^2-r0^2-0.5*a^2)/(a*r0*Ne);
    Iell=I0+Is*sin(theta)+Ic*cos(theta);
    
    % Match info on BCs
    
    solk0.Riphi=ell.RR;
    solk0.Ziphi=ell.ZZ;
    solk0.Iphi=Iell;

    clear source point
    Points =meshData_loc.n;
    source = [ell.RR ell.ZZ] ;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        
        Psi_n_pla=fun_Green_Flux_Loop(source,Points,Iell);
        Psi_n_pla = Psi_n_pla/2/pi;
        
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        
        R_points = Points(:,1);
        Z_points = Points(:,2);
        npt_point = length(R_points);
        
        OPT_PARALLEL = 1;
        n_threads = 24;
        
        R_source = source(:,1);
        Z_source = source(:,2);
        npt_source = length(R_source);
        I_source = Iell;
        
        psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
        Psi_n_pla = psi/2/pi;
        
    end

    
    % %     source.R=ell.RR;
    % %     source.Z=ell.ZZ;
    % %     source.current=ones(size(source.R));
    % %     source.type='ring';
    % %     [GG_p_n]=fun_Field_Loop(source,point);
    % %     GG_p_n.psi = GG_p_n.psi/2/pi;
    % %
    % %     Psi_n_pla = GG_p_n.psi*(Iell);
    
    % First solution
    Psi_n = solk0.Psi_MS + Psi_n_pla;% + 0.3*plaparameter.Ipla);%Iell;
    
% %         figure
% %         edges=pdemesh(meshData.n',meshData.e',[]);
% %         set(edges,'color','k')
% %         set(edges,'linewidth',1)
% %         rr=meshData_loc.n(:,1);
% %         zz=meshData_loc.n(:,2);
% %         axis equal; hold on;
% %         axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
% %         plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k')
% %         rgrid=linspace(Rplot(1),Rplot(2),200);
% %         zgrid=linspace(Zplot(1),Zplot(2),200);
% %         [RR,ZZ] = meshgrid(rgrid,zgrid);
% %         PSI = griddata(rr,zz,Psi_n,RR,ZZ);
% %         contour(RR,ZZ,PSI,100)
    
    
    
    
    %%%%%%%%%
    
elseif SETTINGS.INITIAL_PLASMA_MODEL == 3 % Circular sample plasma with gaussian current density
    
% %     error('SETTINGS.INITIAL_PLASMA_MODEL == 3 no longer available')
    r0=plaparameter.Centroid(1);
    z0=plaparameter.Centroid(2);
    % %     plot(meshData.n(meshData.vess,1),meshData.n(meshData.vess,2),'*')
    
    dist=sqrt((meshData_loc.n(meshData_loc.ind_n_FW,1)-r0).^2+...
        (meshData_loc.n(meshData_loc.ind_n_FW,2)-z0).^2);
%     r_circ=.95*min(dist);
    r_circ=.7*min(dist);
    theta=linspace(0,2*pi,100)';
    circ=[r0+r_circ*cos(theta) z0+r_circ*sin(theta)];
    
% %     plot(circ(:,1),circ(:,2),'*r')
    
    nodes_in = meshData_loc.c_t(find(inpolygon(meshData_loc.c_t(:,1),meshData_loc.c_t(:,2),circ(:,1),circ(:,2))),:);
    
    RR = nodes_in(:,1)-r0;
    ZZ = nodes_in(:,2)-z0;
    beta=pi./r_circ;
    Ipla=1*(cos(beta*sqrt(RR.^2+ZZ.^2))+1);
    
    Ipla = Ipla/sum(Ipla)*plaparameter.Ipla;
    
    source = nodes_in;
    Iell = Ipla;
    Points =meshData_loc.n;

    if SETTINGS.RUN_MEX_ROUTINE == false
        
        Psi_n_pla=fun_Green_Flux_Loop(source,Points,Iell);
        Psi_n_pla = Psi_n_pla/2/pi;
        
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        
        R_points = Points(:,1);
        Z_points = Points(:,2);
        npt_point = length(R_points);
        
        OPT_PARALLEL = 1;
        n_threads = 24;
        
        R_source = source(:,1);
        Z_source = source(:,2);
        npt_source = length(R_source);
        I_source = Iell;
        
        psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
        Psi_n_pla = psi/2/pi;
        
    end

    
    % First solution
    Psi_n = solk0.Psi_MS + Psi_n_pla;% + 0.3*plaparameter.Ipla);%Iell;    
    
    % %         figure
    % %         edges=pdemesh(meshData.n',meshData.e_interface',[]);
    % %         set(edges,'color','k')
    % %         set(edges,'linewidth',1)
    % %         rr=meshData_loc.n(:,1);
    % %         zz=meshData_loc.n(:,2);
    % %         axis equal;
    % %         axis([Rplot Zplot]),  hold on, xlabel('r [m]'), ylabel('z [m]');
    % %         rgrid=linspace(Rplot(1),Rplot(2),200);
    % %         zgrid=linspace(Zplot(1),Zplot(2),200);
    % %         [RR,ZZ] = meshgrid(rgrid,zgrid);
    % %         PSI = griddata(rr,zz,Psi_n,RR,ZZ);
    % %         contour(RR,ZZ,PSI,50)
    
    
end

%%
solk0.Psi=Psi_n;
solk0.Psi_BC=Psi_n(meshData_loc.ind_B);
solk0.Psi_BC_pla = Psi_n_pla(meshData_loc.ind_n_bc);
solk0.Riphi=r0;
solk0.Ziphi=z0;
solk0.Iphi =plaparameter.Ipla;
solk0.r0=r0;
solk0.z0=z0;


