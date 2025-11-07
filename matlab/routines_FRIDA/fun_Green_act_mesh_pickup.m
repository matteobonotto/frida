function [Green_act_pickup] = fun_Green_act_mesh_pickup(meshData,ind_conductors,B_sens,OPT_MEX,OPT_PARALLEL)


BSENS_R = B_sens(:,1);
BSENS_Z = B_sens(:,2);
BSENS_T = B_sens(:,[3 4]);

n_cond = length(ind_conductors);
Green_act_br = zeros(length(BSENS_Z),n_cond);
Green_act_bz = zeros(length(BSENS_Z),n_cond);

Area_all = fun_AreaCoils_FEM(meshData,ind_conductors);

% % tic
for ii=1:n_cond
% %     fprintf('%i\n',ii)
    
% %     ind_t_act_ii = find(meshData.type == ii);
    tri_ii = meshData.t(meshData.type == ind_conductors(ii),:);
% %     nodes_act = meshData.n(tri_ii,:);
    
    filament = [BSENS_R BSENS_Z];
    % define circuits

    nodes_circuit_1 = [meshData.n(tri_ii(:,1),:) ...
        meshData.n(tri_ii(:,2),:) ...
        meshData.n(tri_ii(:,3),:)];
    
    P1 = nodes_circuit_1(:,[1 2]);
    P2 = nodes_circuit_1(:,[3 4]);
    P3 = nodes_circuit_1(:,[5 6]);
    
    % source
    GAUSS_QUAD_DEGREE_circ_1 = 4;
    
    if OPT_MEX == 'F'
        [w_Gauss_circuit_1,P_Gauss_circuit_1,~] = ...
            fun_Gauss_points_triangle_stable(P1,P2,P3,GAUSS_QUAD_DEGREE_circ_1);
    elseif OPT_MEX == 'T'
        [w_Gauss_circuit_1,P_Gauss_circuit_1,~] = ...
            fun_Gauss_points_triangle_stable_mex(P1,P2,P3,GAUSS_QUAD_DEGREE_circ_1);
    end
    
    Area_ii = Area_all(ii);
    
    
    R_source = P_Gauss_circuit_1(:,1);
    Z_source = P_Gauss_circuit_1(:,2);
    npt_source = length(R_source);
    I_source = w_Gauss_circuit_1/Area_ii;
    R_points = filament(:,1);
    Z_points = filament(:,2);
    npt_point = length(R_points);
    
    n_threads = 24;
    
    if OPT_PARALLEL == 'F'
        [out_br,out_bz] = fun_Green_filament_BrBz_SP_f90(npt_source,R_source,Z_source, ...
            I_source,npt_point,R_points,Z_points,0,n_threads);
    elseif OPT_PARALLEL == 'T'
        [out_br,out_bz] = fun_Green_filament_BrBz_SP_f90(npt_source,R_source,Z_source, ...
            I_source,npt_point,R_points,Z_points,1,n_threads);
    end
    
    Green_act_br(:,ii) = out_br;
    Green_act_bz(:,ii) = out_bz;
    
end
% % toc

%%%
vers_r_act = repmat(BSENS_T(:,1),1,n_cond);
vers_z_act = repmat(BSENS_T(:,2),1,n_cond);
Green_act_pickup = Green_act_br.*vers_r_act + Green_act_bz.*vers_z_act;


end


%%
function  [Area_coils]=fun_AreaCoils_FEM(meshData,ind_conductors)

% % figure; hold on
n_cond = length(ind_conductors);
Area_coils = zeros(n_cond,1);
for ii = 1:n_cond
    
    tri_ii = meshData.t(meshData.type == ind_conductors(ii),:);
    
% %     triplot(meshData.t(meshData.type == ii,1:3),meshData.n(:,1),meshData.n(:,2))
    
    P1 = meshData.n(tri_ii(:,1),:);
    P2 = meshData.n(tri_ii(:,2),:);
    P3 = meshData.n(tri_ii(:,3),:);
    
    edge_1 = P2 - P1;
    edge_2 = P3 - P2;

    Area_ii = .5*abs(edge_1(:,1).*edge_2(:,2) - edge_1(:,2).*edge_2(:,1));
    
    Area_coils(ii) = sum(Area_ii);
% %     pause
end

end



