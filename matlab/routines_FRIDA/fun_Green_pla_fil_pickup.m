function Green_pla_pickup = fun_Green_pla_fil_pickup(meshData_loc,B_sens,SETTINGS, OPT_PARALLEL)

%%%
npt_source = size(meshData_loc.n,1);
vers_r = repmat(B_sens(:,3),1,npt_source);
vers_z = repmat(B_sens(:,4),1,npt_source);


%%%
if SETTINGS.RUN_MEX_ROUTINE == false
    
    % %         point.RR = Psi_sens(:,1);
    % %         point.ZZ = Psi_sens(:,2);
    % %         source.R = meshData_loc.n(:,1);
    % %         source.Z = meshData_loc.n(:,2);
    % %         source.current=ones(numel(source.R),1);
    % %         Green_pla_flux_loops2 = fun_GreenMat_Loop_flux(source, point);
    
    I_source = 1;
    
    npt_point = size(B_sens,1);
    
    Green_pla_Br = zeros(npt_point,npt_source);
    Green_pla_Bz = zeros(npt_point,npt_source);
    for ii = 1:npt_source
        
        RZ_source = meshData_loc.n(ii,:);
        [Br_ii,Bz_ii] = fun_Green_BrBz_Loop(RZ_source,B_sens,I_source);
        
        Green_pla_Br(:,ii) = Br_ii;
        Green_pla_Bz(:,ii) = Bz_ii;
        
    end
    
    
else
    
    I_source = 1;
    
    R_points = B_sens(:,1);
    Z_points = B_sens(:,2);
    npt_point = length(R_points);
    
    Green_pla_Br = zeros(npt_point,npt_source);
    Green_pla_Bz = zeros(npt_point,npt_source);
    
    if npt_point < 200
        OPT_PARALLEL = 0;
    end
    
    for ii = 1:npt_source
        
        R_source = meshData_loc.n(ii,1);
        Z_source = meshData_loc.n(ii,2);
        
        [Br_ii,Bz_ii] = fun_Green_filament_BrBz_SP_f90(1, ...
            R_source,...
            Z_source,...
            I_source,...
            npt_point,...
            R_points,...
            Z_points,...
            OPT_PARALLEL,...
            24);
        
        Green_pla_Br(:,ii) = Br_ii;
        Green_pla_Bz(:,ii) = Bz_ii;
        
    end
    
end

%%%
Green_pla_pickup = Green_pla_Br.*vers_r + Green_pla_Bz.*vers_z;






