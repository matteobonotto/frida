function [Green_pas_pickup] = fun_Green_pas_quad_pickup(R_source,Z_source,B_sens,n_Gauss_quad,OPT_PARALLEL)

BSENS_R = B_sens(:,1);
BSENS_Z = B_sens(:,2);
BSENS_T = B_sens(:,[3 4]);

%%%
Green_pas_pickup = zeros(length(BSENS_R),size(R_source,2));
for ii=1:size(R_source,2)
    
    quad = [R_source(:,ii) Z_source(:,ii)];
    [quad_sort,~] = fun_ordinapunti(quad);
    area_quad = polyarea(quad_sort(:,1),quad_sort(:,2));
    j_source_ii = 1./area_quad;
    [out_Br,out_Bz] = fun_Green_quad_BrBz(quad(:,1),quad(:,2),j_source_ii,BSENS_R,BSENS_Z,n_Gauss_quad,OPT_PARALLEL);
    
    out_pickup = sum([out_Br out_Bz].*BSENS_T,2);
    Green_pas_pickup(:,ii) = out_pickup;
    
end