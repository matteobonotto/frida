function [solk0] = fun_Vacuum_FRIDA(meshData,meshData_loc,Conductors,SETTINGS)


n_Gauss = meshData_loc.n_Gauss;
KK_nobc = meshData_loc.KK_nobc;

[Area_coils]=fun_AreaCoils_FEM(meshData,Conductors);


%% Vacuum solution (Integral sources)
if SETTINGS.VAC_METHOD == 1 %
    disp('Vacuum solution (Integral sources) ...')
    
    hh = waitbar(0, 'computing Vacuum solution ...');
    
    Points =meshData_loc.n;
    Points_bc =meshData_loc.n(meshData_loc.ind_n_bc,:);
    
    Psi_MS = zeros(size(Points,1),1);
    Psi_MS_bc = zeros(size(Points_bc,1),1);
    
    for ii=1:Conductors.Nconductors
        waitbar(ii/Conductors.Nconductors)
        
        ind_t_ii = find(meshData.type == ii);
        
        Ic=Conductors.Nturns(ii)*(Conductors.Currents(ii));
        
        Area_coil_ii = Area_coils(ii);
        J_cond = Ic/Area_coil_ii;
            
        ind_t_jj = ind_t_ii(:);
        tri_jj = meshData.t(ind_t_jj,:);
        
        nodes_jj = [meshData.n(tri_jj(:,1),:) ...
            meshData.n(tri_jj(:,2),:) ...
            meshData.n(tri_jj(:,3),:)];
        
        P1 = nodes_jj(:,[1 2]);
        P2 = nodes_jj(:,[3 4]);
        P3 = nodes_jj(:,[5 6]);
        
        if SETTINGS.RUN_MEX_ROUTINE == false
            [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,SETTINGS.GAUSS_QUAD_DEGREE_PLA);
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable_mex(P1,P2,P3,SETTINGS.GAUSS_QUAD_DEGREE_PLA);
        end
        
        j_source = J_cond*ww;
        
        % Psi MS
        clear source point
        point = Points;
        source = nodes;

        if SETTINGS.RUN_MEX_ROUTINE == false
                      
            Psi_MS_bc_ii = fun_Green_Flux_Loop( source, point, j_source);
            Psi_MS_bc_ii = Psi_MS_bc_ii/2/pi;
            
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            
            R_points = point(:,1);
            Z_points = point(:,2);
            npt_point = length(R_points);
            
            OPT_PARALLEL = 1;
            n_threads = 24;
            
            R_source = source(:,1);
            Z_source = source(:,2);
            npt_source = length(R_source);
            I_source = j_source;
            
            psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
            Psi_MS_ii = psi/2/pi;

        end
   
        Psi_MS = Psi_MS + Psi_MS_ii;
        
        
        % Psi MS boundary conditions
        clear source point
        point = Points_bc;
        source = nodes;

        if SETTINGS.RUN_MEX_ROUTINE == false
                      
            Psi_MS_bc_ii = fun_Green_Flux_Loop( source, point, j_source);
            Psi_MS_bc_ii = Psi_MS_bc_ii/2/pi;
            
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            
            R_points = point(:,1);
            Z_points = point(:,2);
            npt_point = length(R_points);
            
            OPT_PARALLEL = 1;
            n_threads = 24;
            
            R_source = source(:,1);
            Z_source = source(:,2);
            npt_source = length(R_source);
            I_source = j_source;
            
            psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
            Psi_MS_bc_ii = psi/2/pi;

        end
   
        Psi_MS_bc = Psi_MS_bc + Psi_MS_bc_ii;
                
        % %         Psi_MS = Psi_MS + Psi_MS_ii;
        % %         Psi_MS_bc_ii=fun_Green_Flux_Loop(nodes,Points_bc,j_source);
        % %         Psi_MS_bc_ii = Psi_MS_bc_ii/2/pi;
        % %
        % %         Psi_MS_bc = Psi_MS_bc + Psi_MS_bc_ii;
        
        % %     figure
        % %     edges=pdemesh(meshData.n',meshData.e',[]);
        % %     set(edges,'color','k')
        % %     set(edges,'linewidth',1)
        % % % %     rr=nodes_flux(:,1); zz=nodes_flux(:,2);
        % %     rr=meshData.n(:,1); zz=meshData.n(:,2);
        % %     axis equal;
        % %     axis([1 4 -1.5 1.5]),  hold on, xlabel('r [m]'), ylabel('z [m]');
        % %     rgrid=linspace(1,4,200);
        % %     zgrid=linspace(-1.5,1.5,200);
        % %     [RR,ZZ] = meshgrid(rgrid,zgrid);
        % %     PSI = griddata(rr,zz,psi_coil_ii,RR,ZZ);
        % %     contour(RR,ZZ,PSI,200);
        % %     PSI = griddata(rr,zz,psi_coil_ii_2,RR,ZZ);
        % %     contour(RR,ZZ,PSI,200);
        % %     colormap('cool');
        % %     title('Vacuum solution')
        
        
    end
    
    close(hh)
    
    meshData_loc.Psi_MS = Psi_MS;
    meshData_loc.Psi_MS_bc = Psi_MS_bc;
    
    solk0.Psi_MS_bc = Psi_MS_bc;
    solk0.Psi_MS = Psi_MS;
    
    
    %% Vacuum solution (FEM + integral BCs)
elseif SETTINGS.VAC_METHOD == 2
    
    disp('Vacuum solution (FEM + integral BCs) ...')
    
    
% %     hh = waitbar(0, 'computing Vacuum solution ...');
    
    Points_bc =meshData_loc.n(meshData_loc.ind_n_bc,:);
    
    Psi_MS_bc = zeros(size(Points_bc,1),1);
    
    for ii=1:Conductors.Nconductors
% %         waitbar(ii/Conductors.Nconductors)
        
        ind_t_ii = find(meshData.type == ii);
        
        Ic=Conductors.Nturns(ii)*(Conductors.Currents(ii));
        
        Area_coil_ii = Area_coils(ii);
        J_cond = Ic/Area_coil_ii;
                
        ind_t_jj = ind_t_ii(:);
        tri_jj = meshData.t(ind_t_jj,:);
        
        nodes_jj = [meshData.n(tri_jj(:,1),:) ...
            meshData.n(tri_jj(:,2),:) ...
            meshData.n(tri_jj(:,3),:)];
        
        P1 = nodes_jj(:,[1 2]);
        P2 = nodes_jj(:,[3 4]);
        P3 = nodes_jj(:,[5 6]);
        
        if SETTINGS.RUN_MEX_ROUTINE == false
            [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,SETTINGS.GAUSS_QUAD_DEGREE_PLA);
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable_mex(P1,P2,P3,SETTINGS.GAUSS_QUAD_DEGREE_PLA);
        end
        
% %         text(sum(nodes(:,1))/numel(nodes(:,1)),sum(nodes(:,2))/numel(nodes(:,1)), ...
% %             num2str(J_cond), 'color', 'black', 'FontSize', 14 );

% %         plot(nodes(:,1),nodes(:,2),'o')
        
        j_source = J_cond*ww;      
        
        clear source point
        point = Points_bc;
        source = nodes;

        if SETTINGS.RUN_MEX_ROUTINE == false
                      
            Psi_MS_bc_ii = fun_Green_Flux_Loop( source, point, j_source);
            Psi_MS_bc_ii = Psi_MS_bc_ii/2/pi;
            
        elseif SETTINGS.RUN_MEX_ROUTINE == true
            
            R_points = point(:,1);
            Z_points = point(:,2);
            npt_point = length(R_points);
            
            OPT_PARALLEL = 1;
            n_threads = 24;
            
            R_source = source(:,1);
            Z_source = source(:,2);
            npt_source = length(R_source);
            I_source = j_source;
            
            psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
            Psi_MS_bc_ii = psi/2/pi;

        end
   
        Psi_MS_bc = Psi_MS_bc + Psi_MS_bc_ii;
        
% %         figure
% %         plot(meshData_loc.n(meshData_loc.ind_n_bc,1),Psi_MS_bc_ii)
% %         pause
        
        % %     figure
        % %     edges=pdemesh(meshData.n',meshData.e',[]);
        % %     set(edges,'color','k')
        % %     set(edges,'linewidth',1)
        % % % %     rr=nodes_flux(:,1); zz=nodes_flux(:,2);
        % %     rr=meshData.n(:,1); zz=meshData.n(:,2);
        % %     axis equal;
        % %     axis([1 4 -1.5 1.5]),  hold on, xlabel('r [m]'), ylabel('z [m]');
        % %     rgrid=linspace(1,4,200);
        % %     zgrid=linspace(-1.5,1.5,200);
        % %     [RR,ZZ] = meshgrid(rgrid,zgrid);
        % %     PSI = griddata(rr,zz,psi_coil_ii,RR,ZZ);
        % %     contour(RR,ZZ,PSI,200);
        % %     PSI = griddata(rr,zz,psi_coil_ii_2,RR,ZZ);
        % %     contour(RR,ZZ,PSI,200);
        % %     colormap('cool');
        % %     title('Vacuum solution')
        
        
    end
% %     close(hh)
    
    if SETTINGS.FIRST_SOLUTION == false
        
        ind_D = meshData_loc.ind_D;
        
        KK_D = KK_nobc(meshData_loc.ind_D,meshData_loc.ind_D);
        KK_bc= KK_nobc(meshData_loc.ind_D,meshData_loc.ind_B);
        
        Sources = zeros(size(ind_D));
        
        Psi_MS_FEM_D = KK_D\(Sources - KK_bc*Psi_MS_bc);
        
        Psi_MS = zeros(meshData_loc.nn,1);
        Psi_MS(meshData_loc.ind_D) = Psi_MS_FEM_D;
        Psi_MS(meshData_loc.ind_B) = Psi_MS_bc;
        
        meshData_loc.Psi_MS = Psi_MS;
        solk0.Psi_MS = Psi_MS;
        
    end
    
    meshData_loc.Psi_MS_bc = Psi_MS_bc;
    solk0.Psi_MS_bc = Psi_MS_bc;
    

    %% Vacuum solution (FEM + integral BCs loaded from INPUT)
elseif SETTINGS.VAC_METHOD == 3 % 
    
    disp('Vacuum solution (FEM + integral BCs loaded from INPUT) ...')
    Psi_MS_bc = SETTINGS.Psi_MS_bc;
    
    ind_D = meshData_loc.ind_D;
    
    KK_D = KK_nobc(meshData_loc.ind_D,meshData_loc.ind_D);
    KK_bc= KK_nobc(meshData_loc.ind_D,meshData_loc.ind_B);
    
    Sources = zeros(size(ind_D));
    
    Psi_MS_FEM_D = KK_D\(Sources - KK_bc*Psi_MS_bc);
    
    Psi_MS = zeros(meshData_loc.nn,1);
    Psi_MS(meshData_loc.ind_D) = Psi_MS_FEM_D;
    Psi_MS(meshData_loc.ind_B) = Psi_MS_bc;
    
    meshData_loc.Psi_MS = Psi_MS;

    solk0.Psi_MS_bc = Psi_MS_bc;
    solk0.Psi_MS = Psi_MS;

end





