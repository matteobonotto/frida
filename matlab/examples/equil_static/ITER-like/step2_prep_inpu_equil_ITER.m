clear
close all
clc
restoredefaultpath


%% Adding folder path
addpath ../../source/routines_FRIDA/


%%
input_Data = './ITER/Data/';

EQUIL_CASE = 1;     % 00080
EQUIL_CASE = 2;     % 00090
EQUIL_CASE = 3;     % 00155
% % EQUIL_CASE = 4;     % 00067
% % EQUIL_CASE = 5;     % 00051
% % EQUIL_CASE = 6;     % 00049

    
input_Equil ='';

% % equilname=090; % ok
equilname=051; % ok


filename_save =  sprintf('_ITER_MR_%1.5i',equilname);  
mesh_Name = 'ITER_meshData_FRIDA_MR';

%% Load equil data
file_psi_data =  [input_Equil  sprintf('Equil_%1.5i_CL4E.mat',equilname)];
load(file_psi_data)

plaparameter.Centroid = globalParameters.currentCentroid;
plaparameter.Ipla     = plasmacurrent;
plaparameter.psibar = fluxFunctionProfiles.normalizedPoloidalFlux;
plaparameter.FdF = fluxFunctionProfiles.ffprime;
plaparameter.dP = fluxFunctionProfiles.staticPprime;


plaparameter.CONFIG = 'TOKAMAK';

plaparameter.equilname = equilname;
plaparameter.input_Equil = input_Equil;

figure(100)
subplot(1,2,1); plot(plaparameter.psibar,plaparameter.FdF,'k-'); hold on;
subplot(1,2,2); plot(plaparameter.psibar,plaparameter.dP,'k-'); hold on

mm = -70:70;
m_tresh = 40;
[psibar,FdF,dP] = fun_ArmonicSmoothing(plaparameter.psibar,plaparameter.FdF,plaparameter.dP,mm,m_tresh);

figure(100)
subplot(1,2,1); plot(psibar,FdF,'r-'); hold on;
subplot(1,2,2); plot(psibar,dP,'r-'); hold on

plaparameter.psibar = psibar;
plaparameter.FdF = FdF;
plaparameter.dP = dP;


%%
plaparameter.CONFIG='TOKAMAK';

%% Save Data
disp(['INPUT_EQUIL' filename_save])


INPUT_EQUIL.meshData=meshData;
INPUT_EQUIL.meshData_loc=meshData_loc;
INPUT_EQUIL.plaparameter=plaparameter;
INPUT_EQUIL.Conductors=Conductors;
INPUT_EQUIL.index=index;
INPUT_EQUIL.index_loc=index_loc;
% % INPUT_INTEGRAL.Psi_MS=Psi_MS;
% % INPUT_INTEGRAL.GG_pla=GG_pla;

% % tic
% % save(['INPUT_Green_' meshName], 'INPUT_PEGASOS')
% % toc


save(['INPUT_EQUIL_FEM' filename_save], 'INPUT_EQUIL')

disp(' ')
disp('Done!!')


return

% % plaparameter.Ipla=plaparameter.Ipla*1.1;

%% Load mesh
% % if N_order == 1
% %     
% %     mesh_name='RFX_meshData_coarse';
% %     % % mesh_name='RFX_meshData_fine';
% %     
% % elseif N_order == 2
% %     
% %     mesh_name = 'RFX_quadratic';
% %     
% % end
N_order = 2;



load([input_Data ,mesh_Name])
keyreg=meshData.type;
nn=size(meshData.n,1);
nt=size(meshData.t,1);
figure
hold on;    grid on
mesh=triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2));
axis equal,  hold on, xlabel('r [m]'), ylabel('z [m]')
set(mesh,'color','k');    
title('RFX mesh with regions')

% Conduttori
ind_conduttori=find(keyreg>0);
    triplot(meshData.t(ind_conduttori,1:3),meshData.n(:,1),meshData.n(:,2),'r')


% Plasma
ind_t_Invess=find(keyreg==-1);
% %     triplot(meshData.t(ind_plasma,:),meshData.n(:,1),meshData.n(:,2),'b')
ind_n_Invess=reshape(meshData.t(ind_t_Invess,:),numel(meshData.t(ind_t_Invess,:)),1);
ind_n_Invess=unique(ind_n_Invess);
pdemesh(meshData.n',[],[meshData.t(ind_t_Invess,1:3)'; ones(1,length(ind_t_Invess))],'EdgeColor','b'); 



%% GEOMETRY
if N_order == 1
 
    P1=meshData.n(meshData.t(:,1),:);  % Primal node 1
    P2=meshData.n(meshData.t(:,2),:);  % Primal node 2
    P3=meshData.n(meshData.t(:,3),:);  % Primal node 3
    % %
    c_t = (P1+P2+P3)/3; % calculate tetrahedron barycenter
    
elseif N_order == 2
    
    P1=meshData.n(meshData.t(:,1),:);  % Primal node 1
    P2=meshData.n(meshData.t(:,2),:);  % Primal node 2
    P3=meshData.n(meshData.t(:,3),:);  % Primal node 3
    P4=meshData.n(meshData.t(:,4),:);  % Primal node 5
    P5=meshData.n(meshData.t(:,5),:);  % Primal node 3
    P6=meshData.n(meshData.t(:,6),:);  % Primal node 3
    % %
    c_t = (P1+P2+P3+P4+P5+P6)/6; % calculate tetrahedron barycenter
    
end

% % plot(c_t(:,1), c_t(:,2), 'ro')

meshData.shape_functions_N_order = N_order;
   
meshData.c_t=c_t;

Rvess=[min(meshData.n(ind_n_Invess,1))...
    max(meshData.n(ind_n_Invess,1))]+[-0.3 0.3];
Zvess=[min(meshData.n(ind_n_Invess,2))...
    max(meshData.n(ind_n_Invess,2))]+[-0.3 0.3];

Rplot=Rvess+[-0.3 0.3];
Zplot=Zvess+[-0.3 0.3];
% % 
meshData.nt=size(meshData.t,1);
meshData.nn=size(meshData.n,1);

meshData.ind_D=setdiff(1:meshData.nn,meshData.b)';
meshData.ind_B=meshData.b;
meshData.nn_D=length(meshData.ind_D);
meshData.nn_B=length(meshData.ind_B);



%% Computational Domain


% % comp_domain = meshData.n(meshData.vess,:);
comp_domain = meshData.n(meshData.ind_n_bc,:);

aa=inpolygon(meshData.c_t(:,1),meshData.c_t(:,2),comp_domain(:,1),comp_domain(:,2));
ind_t_comp=find(aa);
% % ind_t_comp=ind_t_Invess;
meshData_loc.ind_t_comp=ind_t_comp;
aa=meshData.t(ind_t_comp,:)';
ind_n_comp=unique(aa(:),'stable');
plot(meshData.n(ind_n_comp,1),meshData.n(ind_n_comp,2),'ro')

meshData_loc.ind_n_comp=ind_n_comp;
figure
hold on;    grid on; box on;
% % mesh=triplot(meshData.t,meshData.n(:,1),meshData.n(:,2));
edges = pdemesh(meshData.n',meshData.e',[]); 
set(edges,'color','k')
axis equal,  hold on, xlabel('r [m]'), ylabel('z [m]')
% % set(mesh,'color','k');    axis([0.8 4.5 -1.5 1.5])
title('RFX mesh with regions')
% % triplot(meshData.t(ind_t_comp,:),meshData.n(:,1),meshData.n(:,2),'r')
pdemesh(meshData.n',[],[meshData.t(ind_t_comp,1:3)'; ones(1,length(ind_t_comp))],'EdgeColor','r'); 



%% MeshData local
node_new=meshData.n(ind_n_comp,:);

if N_order == 1
    tri_new=zeros(length(ind_t_comp),3);
elseif N_order == 2
    tri_new=zeros(length(ind_t_comp),6);
end

tri_chosen=meshData.t(ind_t_comp,:);
tic
for ii=1:length(ind_n_comp)
    [a,b]=find(tri_chosen==ind_n_comp(ii));
    for jj=1:length(a)
        tri_new(a(jj),b(jj))=ii;
    end
end
toc
pdemesh(node_new',[],[tri_new(:,1:3)'; ones(1,size(tri_new,1))],'EdgeColor','c'); 

meshData_loc.t=tri_new;
meshData_loc.nt = size(meshData_loc.t,1);
meshData_loc.n=node_new;
meshData_loc.type=meshData.type(ind_t_comp);
meshData_loc.P1=[meshData_loc.n(meshData_loc.t(:,1),1), ...
    meshData_loc.n(meshData_loc.t(:,1),2)];
meshData_loc.P2=[meshData_loc.n(meshData_loc.t(:,2),1), ...
    meshData_loc.n(meshData_loc.t(:,2),2)];
meshData_loc.P3=[meshData_loc.n(meshData_loc.t(:,3),1), ...
    meshData_loc.n(meshData_loc.t(:,3),2)];
meshData_loc.c_t=[sum(meshData_loc.P1(:,1)+meshData_loc.P2(:,1)+meshData_loc.P3(:,1),2)/3, ...
    sum(meshData_loc.P1(:,2)+meshData_loc.P2(:,2)+meshData_loc.P3(:,2),2)/3];
meshData_loc.nn=size(meshData_loc.n,1);

ind_t_loc_Invess=find(meshData_loc.type==-1)';
ind_n_loc_Invess=unique(meshData_loc.t(ind_t_loc_Invess,:),'stable');

% % meshData_loc.vess=boundary(meshData_loc.n(ind_n_loc_Invess,1),...
% %     meshData_loc.n(ind_n_loc_Invess,2),.00005);
% % meshData_loc.vess=ind_n_loc_Invess(meshData_loc.vess(1:end-1));

n_points_1 = size(comp_domain,1);

ind_n_loc_Onvess = zeros(n_points_1,1);
for ii = 1:n_points_1
    point_ii = comp_domain(ii,:);
    
    distance = sqrt((point_ii(1)-meshData_loc.n(:,1)).^2 + ...
       (point_ii(2)-meshData_loc.n(:,2)).^2);
    [value,ind_value] = min(distance);
    
% %     if value == 0
        ind_n_loc_Onvess(ii) = ind_value;
% %     end
           
end

meshData_loc.vess = ind_n_loc_Onvess;
% % plot(meshData_loc.n(meshData_loc.vess,1),meshData_loc.n(meshData_loc.vess,2),'*m')



% Area of primal and dual elements (local mesh)
meshData_loc.shape_functions_N_order = N_order;
[Area_tri_loc]=fun_AreaMesh_FEM(meshData_loc);
meshData_loc.Area_tri=Area_tri_loc;

ind_all2loc=zeros(meshData_loc.nn,1);

for ii=1:meshData_loc.nn
   P_loc_ii=meshData_loc.n(ii,:);
   dist_ii=sqrt((P_loc_ii(1)-meshData.n(:,1)).^2 + ...
       (P_loc_ii(2)-meshData.n(:,2)).^2);
   [aa,bb]=min(dist_ii);
   if aa == 0
       ind_all2loc(ii)=bb;
   else
       error('wrong point found')
   end 
   
end

meshData.ind_all2loc=ind_all2loc;
meshData_loc.ind_all2loc=ind_all2loc;

plot(meshData.n(ind_all2loc,1),meshData.n(ind_all2loc,2),'g*')

%% Plasma Domain and boundary conditions
nn=meshData_loc.nn;

ind_FW = zeros(size(meshData.ind_n_FW,1),1);
for ii = 1:length(ind_FW)
    
    point_ii = meshData.n(meshData.ind_n_FW(ii),:);
    distance = sqrt((point_ii(1)-meshData_loc.n(:,1)).^2 + (point_ii(2)-meshData_loc.n(:,2)).^2);
    [aa,bb] = min(distance);
    ind_FW(ii) = bb;
    
end

meshData_loc.ind_n_FW = ind_FW;

ind_n_bc = zeros(size(comp_domain,1),1);
for ii = 1:size(comp_domain,1)
    
    point_ii = comp_domain(ii,:);
    distance = sqrt((point_ii(1)-meshData_loc.n(:,1)).^2 + (point_ii(2)-meshData_loc.n(:,2)).^2);
    [aa,bb] = min(distance);
    ind_n_bc(ii) = bb;
    
end

meshData_loc.ind_n_bc = ind_n_bc;


% % plot(meshData_loc.n(meshData_loc.vess,1),meshData_loc.n(meshData_loc.vess,2),'*m')
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'or')
plot(meshData_loc.n(meshData_loc.ind_n_bc,1),meshData_loc.n(meshData_loc.ind_n_bc,2),'ob')

ind_B=meshData_loc.ind_n_bc;
% % plot(meshData.n(ind_n_comp(ind_B),1),meshData.n(ind_n_comp(ind_B),2),'*b')
% % ind_B=meshData.b;
ind_D=setdiff(1:nn,ind_B);
nn_D=length(ind_D);
nn_B=length(ind_B);

meshData_loc.ind_D=ind_D;
meshData_loc.ind_B=ind_B;
meshData_loc.nn_D=nn_D;
meshData_loc.nn_B=nn_B;


%% Area of triangular elements and nodes
[Area_t]=fun_AreaMesh_FEM(meshData);
meshData.Area_t = Area_t;
% % meshData.Area_n = Area_n;

% % [Area_t,Area_n]=fun_AreaMesh_FEM(meshData_loc);
% % meshData_loc.Area_t = Area_t;
% % meshData_loc.Area_n = Area_n;

%% Conduttori
point.RR=meshData.n(:,1);
point.ZZ=meshData.n(:,2);

Ncond = 12;


% % 
% % Ncond=0;
% % % OH + FS
% % for j=1:Ncircuiti
% %     flag_circuiti=find(KONNAX(j,:));
% %     Ncond=Ncond+numel(flag_circuiti);
% % end
% % % Saddle
% % Ncond=Ncond+numel(find(KONNAX(i_saddle_KONN,:)));

%
Conductors.Nconductors=Ncond;
Conductors.Nturns=zeros(Ncond,1);
Conductors.Currents=zeros(Ncond,1);
Conductors.R=zeros(Ncond,1);
Conductors.Z=zeros(Ncond,1);
Conductors.DR=zeros(Ncond,1);
Conductors.DZ=zeros(Ncond,1);
Conductors.type=zeros(nt,1);


for j=1:Ncond
    
    c=rand(3,1);
    n_turns=coils_data.turns(j);
    Conductors.Currents(j)=coils_data.current(j);
    Conductors.Nturns(j)=coils_data.turns(j);
    ind_t_conductor_j=find(keyreg==j);
    ind_t_conductor{j} = ind_t_conductor_j;
    
    ind_n_coil=unique(reshape(meshData.t(ind_t_conductor_j,:),...
            numel(meshData.t(ind_t_conductor_j,:)),1));
    ind_n_conductor{j}=ind_n_coil;
    pdemesh(meshData.n',[],[meshData.t(ind_t_conductor_j,1:3)'; ones(1,length(ind_t_conductor_j))],'EdgeColor','r'); axis([0.8 4.5 -1.5 1.5]);

% %     plot(meshData.n(ind_n_coil,1),meshData.n(ind_n_coil,2),'ro')
    
    Conductors.R(j)=sum(coils_data.R(j,:))/numel(coils_data.R(j,:));
    Conductors.Z(j)=sum(coils_data.Z(j,:))/numel(coils_data.Z(j,:));
    Conductors.DR(j)=coils_data.DR(j);
    Conductors.DZ(j)=coils_data.DZ(j);

 
end

% % figure
% % for kk=1:Conductors.Nconductors
% %     plot(meshData.n(ind_n_conductor{kk},1),meshData.n(ind_n_conductor{kk},2),'r*');
% %     hold on
% %     plot(meshData.n_D(ind_n_conductor_D{kk},1),meshData.n_D(ind_n_conductor_D{kk},2),'og')
% %     pause
% % end

%%
FW = meshData_loc.n(meshData_loc.ind_n_FW,:);
FW_fill = [FW; FW(1,:)];
comp_domain_fill = [comp_domain; comp_domain(1,:)];

figure
hold on;    grid on; box on;
axis equal,  hold on, 
fill(FW_fill(:,1),FW_fill(:,2),[.07 .62 1],'Linestyle','none')
fill([comp_domain_fill(:,1); FW_fill(:,1)],[comp_domain_fill(:,2); FW_fill(:,2)],'g','Linestyle','none')

tmp = find(meshData.type>0);
sol_tmp = zeros(size(meshData.n(:,2)));

JP.faces=meshData.t(tmp,1:3);
JP.vertices=meshData.n;
JP.facevertexcdata = sol_tmp;

hh=patch(JP,'facecolor','k','edgecolor','none');

plot([comp_domain(:,1); comp_domain(1,1)],[comp_domain(:,2); comp_domain(1,2)],'r','LineWidth',2)


alpha(.6)
fontsize = 10;
set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontsize);
xlabel('r [m]', 'interpreter', 'latex','fontsize',fontsize);
ylabel('z [m]', 'interpreter', 'latex','fontsize',fontsize);

legend('$\Omega_p$','$\Omega_v$','$\Omega_c$','$\partial\Omega$','interpreter', 'latex','fontsize',fontsize-1);

axis([0 13 -8 8])

% % savefig('ITER_with_regions')
% % print('ITER_with_regions', '-depsc', '-r300')

subplot(1,2,1); title('RFX-MOD', 'interpreter', 'latex','fontsize',fontsize);
subplot(1,2,2); title('ITER-like', 'interpreter', 'latex','fontsize',fontsize);

% % savefig('RFX_ITER_with_regions')
% % print('RFX_ITER_with_regions', '-depsc', '-r500')

%% Boundary conditions


Psi_CREATE = load(file_psi_data);

meshData_CREATE.t = t(1:3,:)';
meshData_CREATE.n = p';
meshData_CREATE.e = e';

Psi_CREATE = meshandequil.psi;

meshData_CREATE.t = meshandequil.t(1:3,:)';
meshData_CREATE.n = meshandequil.p';
meshData_CREATE.e = meshandequil.e';

Psi_CREATE = meshandequil.psi;
Psi_CREATE_boundary = meshandequil.psb;

% % Psi_CREATE=2*pi*Psi_CREATE;
disp(' ')
warning('Using poloidal flux per radiant and not the total poloidal flux')
disp(' ')

rr = linspace(3,10,100);
zz = linspace(-6,6,100);
[RR,ZZ] = meshgrid(rr,zz);


figure
      edges=pdemesh(meshData_CREATE.n',meshData_CREATE.e',[]);
      set(edges,'color','k'); set(edges,'linewidth',1); axis equal;
      axis([Rplot Zplot]),  hold on,
      PSI_2 = griddata(meshData_CREATE.n(:,1),meshData_CREATE.n(:,2),Psi_CREATE,RR,ZZ);
      contour(RR,ZZ,PSI_2,100);
% %       contour(RR,ZZ,PSI_2,[-.1178 -.1178],'r')
% %       contour(RR,ZZ,PSI_2,[-.22485 -.22485],'r')
      colormap('cool'); colorbar vert; title('CREATE_L')  
      contour(RR,ZZ,PSI_2,[Psi_CREATE_boundary Psi_CREATE_boundary],'r','LineWidth',2);

      try
% %           switch EQUIL_CASE
% %               case 1
% %                   namefig = 'fig_00080_CREATE';
% %               case 3
% %                   namefig = 'fig_39136_CREATE';
% %               case 4
% %                   namefig = 'fig_00067_CREATE';
% %           end
        namefig=sprintf('fig_%1.5i_CREATE',equilname);
        savefig(namefig)
          
      catch
          
      end


%% Index of triangles
index.ind_t_Invess=ind_t_Invess;
index.ind_n_Invess=ind_n_Invess;
index.ind_n_FW=meshData.ind_n_FW;
% % index.ind_n_VV_int=meshData.ind_n_VV_int;
% % index.ind_n_TSS_int=meshData.ind_n_TSS_int;
index.ind_n_bc=meshData.ind_n_bc;
index.vess=meshData.vess;
index.ind_t_conductor=ind_t_conductor;
index.ind_n_conductor=ind_n_conductor;
% % index.ind_n_conductor_D=ind_n_conductor_D;

% local indexing
index_loc.ind_t_Invess=ind_t_loc_Invess;
index_loc.ind_n_Invess=ind_n_loc_Invess;
index_loc.vess=meshData_loc.vess;
index_loc.ind_n_FW=meshData_loc.ind_n_FW;
index_loc.ind_n_bc=meshData_loc.ind_n_bc;

%% Area coils
[Area_coils]=fun_AreaCoils_FEM(meshData,Conductors,index);
meshData.Area_coils = Area_coils;


%%
plaparameter.CONFIG='TOKAMAK';

%% Save Data
disp(['INPUT_EQUIL' filename_save])


INPUT_EQUIL.meshData=meshData;
INPUT_EQUIL.meshData_loc=meshData_loc;
INPUT_EQUIL.plaparameter=plaparameter;
INPUT_EQUIL.Conductors=Conductors;
INPUT_EQUIL.index=index;
INPUT_EQUIL.index_loc=index_loc;
% % INPUT_INTEGRAL.Psi_MS=Psi_MS;
% % INPUT_INTEGRAL.GG_pla=GG_pla;

% % tic
% % save(['INPUT_Green_' meshName], 'INPUT_PEGASOS')
% % toc


save(['INPUT_EQUIL_FEM' filename_save], 'INPUT_EQUIL')

disp(' ')
disp('Done!!')










