
% % psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
% %
% % % Inputs
% % npt_source = number of source points
% % R_source = vector (npt_source,1)
% % Z_source = vector (npt_source,1)
% % I_source = vector (npt_source,1)
% %
% % npt_point = number of points where psi is computed
% % R_points = vector (npt_point,1)
% % Z_points = vector (npt_point,1)
% %
% % OPT_PARALLEL = 0 for serial execution, 1 for parallel
% % n_threads = number of threads for Open MP
% %
% % % Outputs
% % psi = vector (npt_point,1)

%   source (sources' geometry) - structure including:
%     - R: radial distance form the axis of the coil's centre [m]
%     - Z: vertical distance from z=0 plane of the coil's centre [m]
%
%   point (evaluation points) - structure including:
%     - RR: array of the radial coordinate of the evaluatin points [m]
%     - ZZ: array of the vertical coordinate of the evaluatin points [m]

function out_psi = fun_Green_filament_flux_SP_f90(npt_source,R_source,Z_source, ...
    I_source,npt_point,R_points,Z_points,par,n_threads)

source = [R_source Z_source];
point = [R_points Z_points];
curr = I_source;

out_psi = fun_Green_Flux_Loop( source, point, curr);

