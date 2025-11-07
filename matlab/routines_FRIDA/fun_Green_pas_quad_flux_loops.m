function [Green_pas_flux_loops] = fun_Green_pas_quad_flux_loops(R_source,Z_source,Psi_sens,n_Gauss_quad,OPT_PARALLEL)

PSISENS_R = Psi_sens(:,1);
PSISENS_Z = Psi_sens(:,2);


%%%
Green_pas_flux_loops = zeros(length(PSISENS_R),size(R_source,2));
for ii=1:size(R_source,2)
    
    quad = [R_source(:,ii) Z_source(:,ii)];
    [quad_sort,~] = fun_ordinapunti(quad);
    area_quad = polyarea(quad_sort(:,1),quad_sort(:,2));
    j_source_ii = 1./area_quad;
    out_psi = fun_Green_quad_flux(quad(:,1),quad(:,2),j_source_ii,PSISENS_R,PSISENS_Z,n_Gauss_quad,OPT_PARALLEL);
    Green_pas_flux_loops(:,ii) = out_psi;
    
end