function [b_target] = fun_CalcBrBz_points(meshData_loc,solk1,meshData,Conductors,SETTINGS)

point_b_target = SETTINGS.point_b_target;

point_b_r = point_b_target(:,1);
point_b_z = point_b_target(:,2);


%%

R_points = point_b_r;
Z_points = point_b_z;
npt_point = length(R_points);

OPT_PARALLEL = 1;
n_threads = 24;

R_source = meshData_loc.n(:,1);
Z_source = meshData_loc.n(:,2);
npt_source = length(R_source);
I_source = solk1.Iphi;

[br,bz] = fun_Green_filament_BrBz_SP_f90(npt_source,...
    R_source, ...
    Z_source, ...
    I_source, ...
    npt_point, ...
    R_points, ...
    Z_points, ...
    OPT_PARALLEL, ...
    n_threads);

Bsens_rz_pla = [br bz];



%%%

[Area_coils]=fun_AreaCoils_FEM(meshData,Conductors);

Points_bc = [point_b_r point_b_z];
Bsens_rz_act = zeros(size(Points_bc,1),2);
for ii=1:Conductors.Nconductors
    
    ind_t_ii = find(meshData.type == ii);
    Ic=Conductors.Nturns(ii)*(Conductors.Currents(ii));
    Area_coil = Area_coils(ii);
    J_cond = Ic/Area_coil;
    ind_t_jj = ind_t_ii(:);
    tri_jj = meshData.t(ind_t_jj,:);
    
    nodes_jj = [meshData.n(tri_jj(:,1),:) ...
        meshData.n(tri_jj(:,2),:) ...
        meshData.n(tri_jj(:,3),:)];
    
    P1 = nodes_jj(:,[1 2]);
    P2 = nodes_jj(:,[3 4]);
    P3 = nodes_jj(:,[5 6]);
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable(P1,P2,P3,5);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        [ww,nodes,n_Gauss] = fun_Gauss_points_triangle_stable_mex(P1,P2,P3,5);
    end
    
    I_source_ii = J_cond*ww;
    
    R_points = point_b_r;
    Z_points = point_b_z;
    npt_point = length(R_points);
    
    OPT_PARALLEL = 1;
    n_threads = 24;
    
    R_source = nodes(:,1);
    Z_source = nodes(:,2);
    npt_source = length(R_source);
    I_source = I_source_ii;
    
    [br,bz] = fun_Green_filament_BrBz_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);

    
    Bsens_rz_act_ii = [br bz];    
    Bsens_rz_act = Bsens_rz_act + Bsens_rz_act_ii;
    
end


b_target = Bsens_rz_pla + Bsens_rz_act;









