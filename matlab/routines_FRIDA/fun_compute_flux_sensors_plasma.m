function [psi_sens] = fun_compute_flux_sensors_plasma(meshData_loc,solk,SETTINGS)

fprintf('Computing Psi on flux loops \n')

%% Sensors

rsens_b = SETTINGS.POSTPROCESSING_PSISENS_R;
zsens_b = SETTINGS.POSTPROCESSING_PSISENS_Z;

ind_sel = find(solk.Iphi);

nodes_pla = meshData_loc.n(ind_sel,:);

RUN_MEX_ROUTINE = SETTINGS.RUN_MEX_ROUTINE;
OPT_PARALLEL = SETTINGS.OPT_PARALLEL;

if RUN_MEX_ROUTINE == false
    
    source = nodes_pla;
    curr = solk.Iphi(ind_sel);   
    
    point = [rsens_b zsens_b];
    
    [ psi_sens ] = fun_Green_Flux_Loop( source, point, curr);

    
elseif RUN_MEX_ROUTINE == true
    
    R_source = nodes_pla(:,1);
    Z_source = nodes_pla(:,2);
    I_source = solk.Iphi(ind_sel);
    npt_source = length(R_source);
    
    R_points = rsens_b;
    Z_points = zsens_b;
    npt_point = length(R_points);
    
    n_threads = 24;
    
    if OPT_PARALLEL == false
        out_psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
            I_source,npt_point,R_points,Z_points,0,n_threads);
    elseif OPT_PARALLEL == true
        out_psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
            I_source,npt_point,R_points,Z_points,1,n_threads);
    end
    
    psi_sens = out_psi;

end






end





