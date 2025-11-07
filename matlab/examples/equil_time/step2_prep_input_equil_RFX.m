%% Load scenario data

% For this example, data has already been prepared. Here follows a
% specification on which format to follow.

load('../../data/equil/data_equil_evol.mat')

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









