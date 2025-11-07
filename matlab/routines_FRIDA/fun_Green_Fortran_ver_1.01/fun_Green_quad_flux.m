function [out_psi] = fun_Green_quad_flux(R_source_quad,Z_source_quad,j_source,R_points,Z_points,n_Gauss_quad,OPT_PARALLEL)


if size(R_source_quad,1) ~= 4
    error('source must be a quadrilateral (4x1 vector)')
end

%%

% Gauss Points on Unit square
[xx_norm,yy_norm,ww_norm] = fun_GaussPoints_2D_MB(n_Gauss_quad);
xi = [-1; 1; 1; -1];
eta = [-1; -1; 1; 1];



if size(R_source_quad,2) == 1 % only 1 coil
    
    
    % Map Gauss points onto the real quadrilateral element
    quad = [R_source_quad Z_source_quad];
    [quad_sort,~] = fun_ordinapunti(quad);
    xe = quad_sort(:,1);
    ye = quad_sort(:,2);
    
    P_Gauss_quad = [sum(.25*repmat(xe',n_Gauss_quad^2,1).*(1+xx_norm*xi').*(1+yy_norm*eta'),2) ...
        sum(.25*repmat(ye',n_Gauss_quad^2,1).*(1+xx_norm*xi').*(1+yy_norm*eta'),2)];
    
    area_quad = polyarea(quad_sort(:,1),quad_sort(:,2));
    
    w_Gauss_quad = ww_norm/sum(ww_norm)*area_quad;
    
    R_source = P_Gauss_quad(:,1);
    Z_source = P_Gauss_quad(:,2);
    npt_source = length(R_source);
    npt_point = length(R_points);
    
    I_source = j_source*w_Gauss_quad;
    
    n_threads = 24;
    
    if OPT_PARALLEL == 'F'
        out_psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
            I_source,npt_point,R_points,Z_points,0,n_threads);
    elseif OPT_PARALLEL == 'T'
        out_psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
            I_source,npt_point,R_points,Z_points,1,n_threads);
    end
    
    
    
else % multiple coils
    
    n_coils = size(R_source_quad,2);
    
    out_psi = zeros(length(R_points),1);
    
    for ii=1:n_coils
        
        % Map Gauss points onto the real quadrilateral element
        quad = [R_source_quad(:,ii) Z_source_quad(:,ii)];
        [quad_sort,~] = fun_ordinapunti(quad);
        xe = quad_sort(:,1);
        ye = quad_sort(:,2);
        
        P_Gauss_quad = [sum(.25*repmat(xe',n_Gauss_quad^2,1).*(1+xx_norm*xi').*(1+yy_norm*eta'),2) ...
            sum(.25*repmat(ye',n_Gauss_quad^2,1).*(1+xx_norm*xi').*(1+yy_norm*eta'),2)];
        
        area_quad = polyarea(quad_sort(:,1),quad_sort(:,2));
        
        w_Gauss_quad = ww_norm/sum(ww_norm)*area_quad;
        
        R_source = P_Gauss_quad(:,1);
        Z_source = P_Gauss_quad(:,2);
        npt_source = length(R_source);
        npt_point = length(R_points);
        
        I_source = j_source(ii)*w_Gauss_quad;
        
        n_threads = 24;
        
        if OPT_PARALLEL == 'F'
            out_psi_ii = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
                I_source,npt_point,R_points,Z_points,0,n_threads);
        elseif OPT_PARALLEL == 'T'
            out_psi_ii = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
                I_source,npt_point,R_points,Z_points,1,n_threads);
        end
        
        out_psi = out_psi + out_psi_ii;
                
    end
    
    
    
end




















