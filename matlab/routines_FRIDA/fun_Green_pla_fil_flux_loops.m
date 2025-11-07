function [Green_pla_flux_loops] = fun_Green_pla_fil_flux_loops(meshData_loc,Psi_sens,SETTINGS,OPT_PARALLEL)

if SETTINGS.RUN_MEX_ROUTINE == false
    
    % %         point.RR = Psi_sens(:,1);
    % %         point.ZZ = Psi_sens(:,2);
    % %         source.R = meshData_loc.n(:,1);
    % %         source.Z = meshData_loc.n(:,2);
    % %         source.current=ones(numel(source.R),1);
    % %         Green_pla_flux_loops2 = fun_GreenMat_Loop_flux(source, point);
    
    npt_source = size(meshData_loc.n,1);
    I_source = 1;
    
    npt_point = size(Psi_sens,1);
    
    Green_pla_flux_loops = zeros(npt_point,npt_source);
    for ii = 1:npt_source
        
        RZ_source = meshData_loc.n(ii,:);
        psi_ii = fun_Green_Flux_Loop(RZ_source,Psi_sens,I_source);
        Green_pla_flux_loops(:,ii) = psi_ii;
        
    end
    
else
    
    npt_source = size(meshData_loc.n,1);
    I_source = 1;
    
    R_points = Psi_sens(:,1);
    Z_points = Psi_sens(:,2);
    npt_point = length(R_points);
    
    Green_pla_flux_loops = zeros(npt_point,npt_source);
    
    if npt_point < 200
        OPT_PARALLEL = 0;
    end
    
    for ii = 1:npt_source
        
        R_source = meshData_loc.n(ii,1);
        Z_source = meshData_loc.n(ii,2);
        
        psi_ii = fun_Green_filament_flux_SP_f90(1, ...
            R_source,...
            Z_source,...
            I_source,...
            npt_point,...
            R_points,...
            Z_points,...
            OPT_PARALLEL,...
            24);
        
        Green_pla_flux_loops(:,ii) = psi_ii;
        
    end
    
end

