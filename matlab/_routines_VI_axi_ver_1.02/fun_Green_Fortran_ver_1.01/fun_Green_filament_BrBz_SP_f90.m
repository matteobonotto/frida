
% % [Br,Bz] = fun_Green_filament_BrBz_SP_f90(npt_source,R_source,Z_source,I_source,npt_point,R_points,Z_points,OPT_PARALLEL,n_threads);
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
% % Br = vector (npt_point,1)
% % Bz = vector (npt_point,1)