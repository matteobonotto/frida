clear; close all; clc;

%% INPUT GEOMETRY
theta=linspace(0,2*pi,100)';
RR=2+cos(theta);
ZZ=0+sin(theta);

GEOMETRY=[RR ZZ];

%% WRITE GMSH .geo FILE
% % filename='GMSH_GEOMETRY';
% % [fileID] = write_GMSH(GEOMETRY,0.1,filename);
% % 
% % %% MESH THE GEOMETRY
% % dos(['gmsh -2 ' filename '.geo']);
% % 
% % %% READ MESH
% % [p,e,t]=readGMSH(filename);
% % 
% % figure; hold on, axis equal;
% % triplot(t,p(:,1),p(:,2));

%% AUTOMATIC MESH WITH GMSH
tic
filename='GMSH_GEOMETRY';
SPACING=.01;
[p,e,t]=autoMESH(filename,GEOMETRY,SPACING);
figure; hold on, axis equal;
triplot(t,p(:,1),p(:,2));
toc
