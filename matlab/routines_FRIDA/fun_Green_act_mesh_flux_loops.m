function [Green_act_flux_loops] = fun_Green_act_mesh_flux_loops(meshData,ind_conductors,Psi_sens,OPT_MEX,OPT_PARALLEL)

n_cond = length(ind_conductors);
Green_act_flux_loops = zeros(size(Psi_sens,1),n_cond);

for ii=1:n_cond
    % %     fprintf('%i\n',ii)
    
    tri_ii = meshData.t(meshData.type == ind_conductors(ii),:);
    
    temp_M = fun_M_Green_Axi_Arbitrary_Mesh_filament(meshData.n,tri_ii,Psi_sens,OPT_MEX,OPT_PARALLEL);
    Green_act_flux_loops(:,ii) = temp_M;
    
end

