

clc; clearvars; close all;


restoredefaultpath
dir_FRIDA = '../../source/routines_FRIDA/';
addpath(genpath(dir_FRIDA))
 
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% geometry per evol
% meshData_ext -> struct containing mesh data of active and passive
%                 conductors. Can be generated via the following script 
%                 .geometry_generation/Matlab2Gmsh_RFX_VI.m
% meshdata_pla -> struct containing mesh data of plasma region only. 
%                 Can be generated via the following script 
%                 .geometry_generation/Matlab2Gmsh_RFX_plasma.m

input_Data = '../../data/';

load([input_Data, 'geo/', 'RFXmod_mesh_VI_linear.mat'])
% % meshData_ext = meshData;
load([input_Data, 'geo/', 'RFXmod_mesh_pla.mat'])



%% Geometry
% Geo info (meshData_pla and meshData_ext) 

% Save meshData_pla and meshData_ext in the same folder where this main.m is located
filename_geo = './data_in_FRIDA/INPUT_FRIDA_geo_RFXmod_36922_rampdown'; % personalize suffix of geo data file and save it
save(filename_geo,'meshData_pla', 'meshData_ext')



%% Magnetic sensors

%%%%% flux-loops
theta_fp =  [22.5 67.5 112.5 157.5 202.5 247.5 292.5 337.5];
% % theta_fp =  linspace(0,360,100); theta_fp = theta_fp(1:end-1);
thetarad = theta_fp*pi/180;
r_fluxloop = 0.5065;
PSISENS_R = [1.995 + r_fluxloop*cos(thetarad)]';
PSISENS_Z = [r_fluxloop*sin(thetarad)]';

flux_loops = [PSISENS_R PSISENS_Z];

SETTINGS.flux_loops = flux_loops;


%%%%% pick-up
theta = [27,72,117,162,207,252,297,342];
% % theta =  linspace(0,360,100); theta = theta(1:end-1);
thetarad=theta*pi/180;
r_pickup = 0.5085;
BSENS_R = [1.995 + r_pickup*cos(thetarad)]';
BSENS_Z = [r_pickup*sin(thetarad)]';
BSENS_T = [-sin(thetarad)' cos(thetarad)'];

pickup = [BSENS_R BSENS_Z BSENS_T];

SETTINGS.pickup = pickup;

% % phi_0 = 70*pi/180;
% % phi = linspace(0,2*pi,25) - phi_0;
% % phi = phi(1:end-1);
% % 
% % XFLD_mat = BSENS_R*cos(phi);
% % YFLD_mat = BSENS_R*sin(phi);
% % ZFLD_mat = repmat(BSENS_Z,1,length(phi));
% % 
% % XFLD = XFLD_mat(:);
% % YFLD = YFLD_mat(:);
% % ZFLD = ZFLD_mat(:);
% % 
% % figure
% % plot_mesh_CARIDDI('RFX_active_passive_renumber.msh',3,[.5 .5 .5])
% % % % plot_mesh_CARIDDI('RFX_active_passive_renumber.msh')
% % hold on
% % plot3(XFLD(1:8),YFLD(1:8),ZFLD(1:8),'bo', 'LineWidth',2)
% % 
% % fprintf('%f, ',XFLD.')
% % 
% % meshData_ext.n(:,1)
% % set(mesh_vac,'color',0*[1 1 1]);


%% Load scenario data

load([input_Data, 'equil/data_equil_evol.mat'])

% ----------------------------------------------- 
% IE_evol (struct containing equil info during discharge)
%     plaparameter: [1×1 struct] -> plasma info
%       Conductors: [1×1 struct] -> active coils' info
%         time_sim: [1×200 double] -> time axis
%         
%         
% IE_evol.plaparameter
%         Ipla: [1×200 double] -> plasma current during discharge
%       beta_0: [1×200 double] -> 1st current density param during discharge
%      alpha_M: [1×200 double] -> 2nd current density param during discharge
%      alpha_N: [1×200 double] -> 3rd current density param during discharge
% 
% [in case we are using ff' and p' profiles rather than [beta_0, alpha_M,
% alpha_N] and defined, at each instant, by a vector of 101 points, we'll have
%      
% IE_evol.plaparameter
%         Ipla: [1×200 double] -> plasma current during discharge
%          FdF: [101×200 double] -> ff' profile during discharge
%           dP: [101×200 double] -> p' profile during discharge
%
%
% IE_evol.Conductors
%     Nconductors: 64 -> total number of conductors
%          Nturns: [64×1 double] -> number of turns per each conductor
%        Currents: [64×200 double] -> coils' currents during discharge
% ----------------------------------------------- 


% Save IE_evol in the same folder where this main.m is located
Conductors = IE_evol.Conductors;
Conductors.time_sim = IE_evol.time_sim;
plaparameter = IE_evol.plaparameter;
plaparameter.time_sim = IE_evol.time_sim;
filename_equil = './data_in_FRIDA/INPUT_FRIDA_equil_RFXmod_36922_rampdown';
save(filename_equil,...
    'plaparameter', ...
    'Conductors')

% For evol sim we need also initial conditions (i.e., at t=0):
%   J_c_0 -> current density on passive conductors. If not available, don't 
%            specify anything and it will be automatically set to zero
%   v_gap_t0 -> voltage between gaps (that may be present if the passive 
%               conductors are toroidally not continuous). Usually not
%               considered for common tokamaks. If not available, don't 
%               specify anything and it will be automatically set to zero



%% SETTINGS
filename_out = 'temp_out';
filename_preproc = 'RFXmod_36922_rampdown'; % specify preproc data filename (if preproc already run)
filename_save = 'temp_save';

SETTINGS.filename_geo = filename_geo;
SETTINGS.filename_equil = filename_equil;
SETTINGS.filename_preproc = filename_preproc;

SETTINGS.filename_save = filename_save;
SETTINGS.filename_out = filename_out;

SETTINGS.FIRST_SOLUTION = false;
SETTINGS.SAVE_OUTPUT = false;
SETTINGS.SOLVER = 'NR';
SETTINGS.RUN_MEX_ROUTINE = true;

SETTINGS.VAC_METHOD = 2;

SETTINGS.plaparameter.Centroid = [2 -.05];
SETTINGS.FIGURES_DEBUG = 0;

% SETTINGS.J_c_t0 = J_c_t0;
% SETTINGS.v_gap_t0 = v_gap_t0;

SETTINGS.QuasiNewton = true;
SETTINGS.QuasiNewton_factor = 0.5;

SETTINGS.ii_START = 1   ;
SETTINGS.PREPROC = true;
SETTINGS.FIGURES = true;

% Run type:
% 0 -> static vacuum
% 1 -> static plasma
% 2 -> evol vacuum
% 3 -> evol plasma
SETTINGS.RUN = 2; 

%% Run FRIDA and store the results
close all
[OUT_FRIDA_TD] = main_run_FRIDA_evol_VI(SETTINGS);

save('OUT_FRIDA_TD_RFXmod_36922_rampdown_11.mat', ...
    'OUT_FRIDA_TD' ,...
    'IE_evol', ...
    'meshData_pla', ...
    'meshData_ext');








