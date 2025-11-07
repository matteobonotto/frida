%% MATLAB 2 GMSH x RFX MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Write RFX geometrical data write in a suitable format for GMSH mesh
%   generation.
%   
%   Matteo Bonotto 05/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear all;  clc;  close all;

%% Input Geometry
input_Data='Data/';
load([input_Data ,'RFX_create_l_saddle_01.mat'],'keyreg');
load([input_Data ,'RFX_create_l_saddle_01_linear.mat']);
load([input_Data ,'kon.mat']);

%%
nt=size(meshData.t,1);
P1=[meshData.n(meshData.t(:,1),1) meshData.n(meshData.t(:,1),2)];
P2=[meshData.n(meshData.t(:,2),1) meshData.n(meshData.t(:,2),2)];
P3=[meshData.n(meshData.t(:,3),1) meshData.n(meshData.t(:,3),2)];
centro_t=(P1+P2+P3)/3;
figure
    hold on;    grid on
% %     pdemesh(p,[],t)
    mesh=triplot(meshData.t,meshData.n(:,1),meshData.n(:,2),'k');
    axis equal,  hold on, xlabel('r [m]'), ylabel('z [m]')
    axis([0.8 4.5 -1.5 1.5])
% Conduttori
ind_conduttori=find(keyreg>0);
%     triplot(meshData.t(ind_conduttori,:),meshData.n(:,1),meshData.n(:,2),'r')
% shell 
ind_shell=find(keyreg>=117 & keyreg<=176);
    triplot(meshData.t(ind_shell,:),meshData.n(:,1),meshData.n(:,2),'g')
% vessel
ind_vessel=find(keyreg>=57 & keyreg<=116);
    triplot(meshData.t(ind_vessel,:),meshData.n(:,1),meshData.n(:,2),'y')
% TSS
ind_TSS=find(keyreg>=177 & keyreg<=236);
    triplot(meshData.t(ind_TSS,:),meshData.n(:,1),meshData.n(:,2),'b')
% %     plot(centro_t(ind_TSS,1),centro_t(ind_TSS,2),'bo')
% saddle
ind_saddle=find(keyreg>=237 & keyreg<=244);
    triplot(meshData.t(ind_saddle,:),meshData.n(:,1),meshData.n(:,2),'c')
% %     plot(centro_t(ind_saddle,1),centro_t(ind_saddle,2),'co')
% Aactive coils (256 conduttori collegati in 194 circuiti)
ind_active=find(keyreg>=1 & keyreg<=56);
    triplot(meshData.t(ind_active,:),meshData.n(:,1),meshData.n(:,2),'r')
%     
%% VV, Shell, TSS 
n_points=250;
VV_in=[1.995+0.475*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
    0.47*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
VV_ex=[1.995+0.505*sin(linspace(0,2*pi-2*pi/n_points,n_points))' ...
    0.5*cos(linspace(0,2*pi-2*pi/n_points,n_points))'];

n_points=200;
Shel_in=[1.995+0.5115*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
    0.506*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
Shel_ex=[1.995+0.5145*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
    0.509*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];

TSS_in=[2+0.553*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
    0.553*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];
TSS_ex=[2+0.6*cos(linspace(0,2*pi-2*pi/n_points,n_points))' ...
    0.6*sin(linspace(0,2*pi-2*pi/n_points,n_points))'];

plot(VV_in(:,1),VV_in(:,2),'b.')
plot(VV_ex(:,1),VV_ex(:,2),'b.')
plot(Shel_in(:,1),Shel_in(:,2),'b.')
plot(Shel_ex(:,1),Shel_ex(:,2),'b.')
plot(TSS_in(:,1),TSS_in(:,2),'b.')
plot(TSS_ex(:,1),TSS_ex(:,2),'b.')

%% Geometria Saddle coils
conduttori={};
% % for ii=1:numel(ind_saddle)
% %     triplot(meshData.t(ind_saddle(ii),:),meshData.n(:,1),meshData.n(:,2),'k')
% %     pause
% % end
ii=1;
for j=1:4
    ind_saddle_coil=ind_saddle((j-1)*4+1:j*4);
    ind_n_saddle=reshape(meshData.t(ind_saddle_coil,:),...
        numel(meshData.t(ind_saddle_coil,:)),1);
    ind_n_saddle=unique(ind_n_saddle);
    coils.rr=meshData.n(ind_n_saddle,1);
    coils.zz=meshData.n(ind_n_saddle,2);
    [coils_sort]=fun_ordinapunti([coils.rr coils.zz]);
    plot(coils.rr,coils.zz,'.r')
% %     pause
    conduttori{ii}=coils_sort;
    ii=ii+1;
end
% conduttori={} ii=1->4 saddle

%% Geometria Attivi
load([input_Data ,'RFX_create_l_saddle_01.mat']); 
meshData.t=t(1:6,:)';
meshData.n=p';
ii=5;
% % for j=1:14
% %     c=rand(3,1);
% %     flag_circuiti=find(KONNAX(j,:));
% %     ind_circuito=[];
% %     for k=1:numel(flag_circuiti)
% %         ind_t_subcirc=find(keyreg==flag_circuiti(k));
% %         ind_circuito=[ind_circuito ind_t_subcirc];
% %         % Green Matrix for subcoil of the circuit
% %         ind_n_subcoil=unique(reshape(meshData.t(ind_t_subcirc,:),...
% %             numel(meshData.t(ind_t_subcirc,:)),1));
% %         r_coil=meshData.n(ind_n_subcoil,1);
% %         z_coil=meshData.n(ind_n_subcoil,2);
% % % %         coils.rr=[min(r_coil); max(r_coil); max(r_coil); min(r_coil)]; 
% % % %         coils.zz=[min(z_coil); min(z_coil); max(z_coil); max(z_coil)]; 
% %         s=0.2;
% %         if j==12
% %             if k==1
% %                 coil_j12_k1=[2.378000000000000   0.669000030000000;...
% %                     2.378000000000000   0.580999970000000;...
% %                     2.415999900000000   0.580999970000000;...
% %                     2.415999900000000   0.541000010000000;...
% %                     2.448000000000000   0.541000010000000;...
% %                     2.448000000000000   0.629000010000000;...
% %                     2.428999900000000   0.629000010000000;...
% %                     2.428999900000000   0.643999990000000;...
% %                     2.410000100000000   0.643999990000000;...
% %                     2.410000100000000   0.653999980000000;...
% %                     2.391000000000000   0.653999980000000;...
% %                     2.391000000000000   0.669000030000000];
% %                 coils.rr=coil_j12_k1(:,1);
% %                 coils.zz=coil_j12_k1(:,2);
% %             elseif k==2
% %                 coil_j12_k2=[2.378000000000000   -0.669000030000000;...
% %                     2.378000000000000   -0.580999970000000;...
% %                     2.415999900000000   -0.580999970000000;...
% %                     2.415999900000000   -0.541000010000000;...
% %                     2.448000000000000   -0.541000010000000;...
% %                     2.448000000000000   -0.629000010000000;...
% %                     2.428999900000000   -0.629000010000000;...
% %                     2.428999900000000   -0.643999990000000;...
% %                     2.410000100000000   -0.643999990000000;...
% %                     2.410000100000000   -0.653999980000000;...
% %                     2.391000000000000   -0.653999980000000;...
% %                     2.391000000000000   -0.669000030000000];
% %                 coils.rr=coil_j12_k2(:,1);
% %                 coils.zz=coil_j12_k2(:,2);
% %             end
% %         elseif j==13
% %             if k==1
% %                 coil_j13_k1=[2.559000000000000   0.470999990000000;...
% %                     2.559000000000000   0.382999990000000;...
% %                     2.596999900000000   0.382999990000000;...
% %                     2.596999900000000   0.342999990000000;...
% %                     2.628999900000000   0.342999990000000;...
% %                     2.628999900000000   0.430999990000000;...
% %                     2.609999900000000   0.430999990000000;...
% %                     2.609999900000000   0.446000010000000;...
% %                     2.591000100000000   0.446000010000000;...
% %                     2.591000100000000   0.456000000000000;...
% %                     2.572000000000000   0.456000000000000;...
% %                     2.572000000000000   0.470999990000000];
% %                 coils.rr=coil_j13_k1(:,1);
% %                 coils.zz=coil_j13_k1(:,2);
% %             elseif k==2
% %                 coil_j13_k2=[2.559000000000000   -0.470999990000000;...
% %                     2.559000000000000   -0.382999990000000;...
% %                     2.596999900000000   -0.382999990000000;...
% %                     2.596999900000000   -0.342999990000000;...
% %                     2.628999900000000   -0.342999990000000;...
% %                     2.628999900000000   -0.430999990000000;...
% %                     2.609999900000000   -0.430999990000000;...
% %                     2.609999900000000   -0.446000010000000;...
% %                     2.591000100000000   -0.446000010000000;...
% %                     2.591000100000000   -0.456000000000000;...
% %                     2.572000000000000   -0.456000000000000;...
% %                     2.572000000000000   -0.470999990000000];
% %                 coils.rr=coil_j13_k2(:,1);
% %                 coils.zz=coil_j13_k2(:,2);
% %             end     
% %         else
% %             ind_bound=boundary(r_coil,z_coil,s);
% %             ind_bound=ind_bound(1:end-1);
% %             coils.rr=r_coil(ind_bound);
% %             coils.zz=z_coil(ind_bound);
% %         end            
% % % %             plot(r_coil,z_coil,'*','color',c)
% % % %             plot(r_coil,z_coil,'.','color',c)
% %         if j==5 || j==6
% %             V1=max(coils.rr);
% %             V2=min(coils.rr);
% %             ind1=find(coils.rr<0.5*(V1+V2));
% %             ind2=find(coils.rr>0.5*(V1+V2));
% %             coils1.rr=coils.rr(ind1);
% %             coils1.zz=coils.zz(ind1);
% %             [coils1_sort]=fun_ordinapunti([coils1.rr coils1.zz]);
% %             coils2.rr=coils.rr(ind2);
% %             coils2.zz=coils.zz(ind2);
% %             [coils2_sort]=fun_ordinapunti([coils2.rr coils2.zz]);
% %             plot(coils1.rr,coils1.zz,'oc')
% %             plot(coils2.rr,coils2.zz,'or')
% %             conduttori{ii}=coils1_sort;
% %             conduttori{ii+1}=coils2_sort;
% %             cond{ii-4}=[min(coils1_sort(:,1)) min(coils1_sort(:,2));...
% %                 max(coils1_sort(:,1)) min(coils1_sort(:,2));...
% %                 max(coils1_sort(:,1)) max(coils1_sort(:,2));...
% %                 min(coils1_sort(:,1)) max(coils1_sort(:,2))];
% %             cond{ii+1-4}=[min(coils2_sort(:,1)) min(coils2_sort(:,2));...
% %                 max(coils2_sort(:,1)) min(coils2_sort(:,2));...
% %                 max(coils2_sort(:,1)) max(coils2_sort(:,2));...
% %                 min(coils2_sort(:,1)) max(coils2_sort(:,2))];
% %             ii=ii+2;
% %         elseif j==12 || j==13
% %             cond{ii-4}=[coils.rr coils.zz];  
% %             ii=ii+1;
% %         else
% %             plot(coils.rr,coils.zz,'*','color',c)
% %             [coils_sort]=fun_ordinapunti([coils.rr coils.zz]);
% %             conduttori{ii}=coils_sort;
% %             cond{ii-4}=[min(coils_sort(:,1)) min(coils_sort(:,2));...
% %                 max(coils_sort(:,1)) min(coils_sort(:,2));...
% %                 max(coils_sort(:,1)) max(coils_sort(:,2));...
% %                 min(coils_sort(:,1)) max(coils_sort(:,2))];
% %             ii=ii+1;
% %         end
% %     end
% % % %     pause
% % end

load([input_Data ,'RFX_Geometria_Conduttori.mat']);
for ii=1:size(cond,2)
% %     conduttore=cond{ii};
% %     plot(conduttore(:,1),conduttore(:,2),'*k')
% %     pause
    conduttori{ii+4}=cond{ii};
end
%%
box=[3 -0.9; 4 -0.9; 4 0.9; 3 0.9];

%% Boundary
boundary=[50*cos(sort(linspace(-0.5*pi,0.5*pi,10),'descend'))' ...
    50*sin(sort(linspace(-0.5*pi,0.5*pi,10),'descend'))']; 
boundary(1,:)=[0 50];
boundary(end,:)=[0 -50];
boundary=[0 1.5; boundary; 0 -1.5];
plot(boundary(:,1),boundary(:,2),'.')
% % axis([-1 51 -50 50]) 

%% Shapes per GMSH
shape={};
shape{1}=VV_in;
shape{2}=VV_ex;
shape{3}=Shel_in;
shape{4}=Shel_ex;
shape{5}=TSS_in;
for ii=1:4
    ii+5;
    shape{ii+5}=conduttori{ii};
end
shape{10}=TSS_ex;
for ii=5:numel(conduttori)
    ii+1+5;
    shape{ii+1+5}=conduttori{ii};
end
shape{numel(shape)+1}=box;
shape{numel(shape)+1}=boundary;

n_point=1;
n_line=1;
n_loop=1;
spacing_factor=ones(size(shape,2),1);
spacing_factor(1)=0.015;       % plasma, vessel (in)
spacing_factor(2)=0.025;       % vessel (out)
spacing_factor(3)=0.015;       % shell (in)
spacing_factor(4)=0.01;        % shell (out)
spacing_factor(5)=0.03;        % TSS (in)
spacing_factor(6:9)=0.02;      % Saddle
spacing_factor(10)=0.05;       % TSS (out)
spacing_factor(11:end-1)=0.02; % Conduttori
spacing_factor(end-1)=0.06;    % Box   
spacing_factor(end)=10;        % Boundary
for ii=1:size(shape,2)-1
    refinement{ii}=spacing_factor(ii)*ones(size(shape{ii},1),1);
end
refinement{ii}=[0.03 0.1 0.1 0.03];
refinement{ii+1}=[0.2; spacing_factor(ii+1)*ones(size(shape{ii+1},1)-2,1); 0.2];

% Points
if exist('GMSHgeometry.geo', 'file')==2
  delete('GMSHgeometry.geo');
end
diary('GMSHgeometry.geo')
for ii=1:size(shape,2)
    element=shape{ii};
    spacing=refinement{ii};
    for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        n_point=n_point+1;
    end
    
    disp(' ')
    
    % Line
    for i=1:size(element,1)
        if i<size(element,1)
            disp(['Line(',num2str(n_line),')={', num2str(n_line),...
                ',', num2str(n_line+1), '};'])
        elseif i==size(element,1)
            disp(['Line(',num2str(n_line),')={', num2str(n_line),...
                ',', num2str(n_line-size(element,1)+1), '};'])
        end
        n_line=n_line+1;
    end
    
    disp(' ')
    
    % Line Loop
    disp(['Line Loop(',num2str(n_loop),')={'])
    for i=1:size(element,1)
        if i<size(element,1)
            disp([num2str(n_line-size(element,1)+i-1), ','])
        elseif i==size(element,1)
            disp(num2str(n_line-size(element,1)+i-1))
        end
    end
    n_loop=n_loop+1;
    disp('};')
    
    disp(' ')

end

disp(' ')
% Plane surface
disp(['Plane Surface(',num2str(1), ')={', num2str(1),'};'])
disp(['Plane Surface(',num2str(2), ')={', num2str(2),', ', num2str(1),'};'])
disp(['Plane Surface(',num2str(3), ')={', num2str(3),', ', num2str(2),'};'])
disp(['Plane Surface(',num2str(4), ')={', num2str(4),', ', num2str(3),'};'])
disp(['Plane Surface(',num2str(5), ')={', num2str(5),', ', num2str(4),'};'])
disp(['Plane Surface(',num2str(6), ')={', num2str(6),'};'])
disp(['Plane Surface(',num2str(7), ')={', num2str(7),'};'])
disp(['Plane Surface(',num2str(8), ')={', num2str(8),'};'])
disp(['Plane Surface(',num2str(9), ')={', num2str(9),'};'])
disp(['Plane Surface(',num2str(10), ')={', num2str(10),', ',...
    num2str(9),', '...
    num2str(8),', '...
    num2str(7),', '...
    num2str(6),', '...
    num2str(5),'};'])
n_surf=11;
ind_surf=11;
for ii=10:size(shape,2)-2
    disp(['Plane Surface(',num2str(n_surf), ')={', num2str(ind_surf),'};'])
    n_surf=n_surf+1;
    ind_surf=ind_surf+1;    
end

disp(['Plane Surface(',num2str(n_surf), ')={', num2str(n_loop-1),', '])
for ii=0:59+1
    disp([num2str(ii+9+1),', '])
end
disp([num2str(ii+9+1+1), '};'])



diary off


