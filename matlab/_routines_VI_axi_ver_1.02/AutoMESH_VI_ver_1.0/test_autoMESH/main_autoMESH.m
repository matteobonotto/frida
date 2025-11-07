close all; clc; clear;

%%
load SAMPLE_GEO

figure; hold on; axis equal
plot(GEOMETRY(:,1),GEOMETRY(:,2),'.-')

cd fun_autoMESH
filename='stab_plate';
SPACING=.1;
tic
disp('MESHING GEOMETRY')
[p,e,t]=autoMESH(filename,GEOMETRY,SPACING);
toc
triplot(t,p(:,1),p(:,2))
cd ../
