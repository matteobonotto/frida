%------------------------------------------------------------------------%
%------ Gmsh to Matlab script: Import mesh to matlab---------------------%
%------------------------------------------------------------------------%

clc
close all
clear 

%-----------------------------------------------------------------------%
% dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range 
% bounded by row offsets R1 and R2 and column offsets C1 and C2.
%-----------------------------------------------------------------------%
file    =  ('RFX_meshData_1.msh');

% no of nodes is mentioned in 5th row and first column

N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

node_id     = dlmread(file,'',[5 0 4+N_n 0]);
nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);

%------- 2D Geometry

two_d_nodes = nodes(:,1:2);
elem_type   = elements(:,2);

%--- find the starting indices of 2D elements
two_ind = 1;
for i = 1:N_e
    if(elem_type(i) ~= 2)
        two_ind = two_ind+1;
    end
end
%----------------------------------------------

two_d_elements(1:N_e-two_ind,1:3) = 0;
k = 1;
for i = two_ind:N_e
   two_d_elements(k,1:3) = elements(i,6:8);
   k = k+1;
end

% edges
% % tri_srt=two_d_elements;
% % ntri=size(two_d_elements,1);
% % lati=zeros(1,2);
% % lati_tot=zeros(3*size(two_d_elements,1),2);
% % for i=1:ntri
% %     triangle=tri_srt(i,:);
% %     lati_loc=[triangle(1) triangle(2);
% %         triangle(2) triangle(3);
% %         triangle(3) triangle(1)];
% %     if i==1
% %         lati=lati(2:end,:);
% %     end
% %     lati_tot(3*i-2:3*i,:)=lati_loc;
% % end
% % lati_tot=sort(lati_tot,2);
% % [lati,~,~]=unique(lati_tot,'rows');
% % 
% % ind_lati_tri=zeros(ntri,3);
% % var_loc=zeros(size(lati));
% % tic
% % parfor i=1:ntri
% %     triangolo=tri_srt(i,:);
% %     var_loc=repmat(sort([triangolo(2) triangolo(1)]),size(lati,1),1);
% %     pos1=find(all(lati==var_loc,2));
% %     var_loc=repmat(sort([triangolo(3) triangolo(2)]),size(lati,1),1);
% %     pos2=find(all(lati==var_loc,2));
% %     var_loc=repmat(sort([triangolo(1) triangolo(3)]),size(lati,1),1);
% %     pos3=find(all(lati==var_loc,2));
% %     ind_lati=[pos1 pos2 pos3];
% %     ind_lati_tri(i,:)=ind_lati;
% % end
% % toc

% % load ITER_meshdata_CREATE; meshDataTOT=meshData;
% % two_d_elements=meshDataTOT.t(1:3,:)';
% % two_d_nodes=meshDataTOT.p';
%---- visualize in matlab ---------------------

figure(1)
c=[0.5 0.5 0.5];
triplot(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('GMsh to MATLAB import','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white'); 
axis equal
axis([0.8 4.5 -1.5 1.5])
hold on

% % plot(two_d_nodes(two_d_elements(1:280,1),1),...
% %     two_d_nodes(two_d_elements(1:280,1),2),'r*');
% % hold on

input_Data='Data/';
load([input_Data ,'RFX_Geometria_Conduttori.mat']);
load([input_Data ,'RFX_create_l_saddle_01_linear.mat']);
load([input_Data ,'kon.mat']);
% % for ii=1:size(cond,2)
% %     conduttore=cond{ii};
% %     plot(conduttore(:,1),conduttore(:,2),'*c')
% % end

%baricentro triangoli
P1=[two_d_nodes(two_d_elements(:,1),1), ...
    two_d_nodes(two_d_elements(:,1),2)];
P2=[two_d_nodes(two_d_elements(:,2),1), ...
    two_d_nodes(two_d_elements(:,2),2)];
P3=[two_d_nodes(two_d_elements(:,3),1), ...
    two_d_nodes(two_d_elements(:,3),2)];
centro=[sum(P1(:,1)+P2(:,1)+P3(:,1),2)/3, ...
    sum(P1(:,2)+P2(:,2)+P3(:,2),2)/3];

% % plot(baricentro(:,1),baricentro(:,2),'.r');

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

type=zeros(size(two_d_elements,1),1);
% elementi interni al vessel
[ind_t_pla,~]=inpolygon(centro(:,1),centro(:,2),...
    VV_in(:,1),VV_in(:,2));
ind_t_pla=find(ind_t_pla);
type(ind_t_pla)=-1;
    triplot(two_d_elements(ind_t_pla,:),two_d_nodes(:,1),two_d_nodes(:,2),'b')
% Vessel    
[ind_t_vess,~]=inpolygon(centro(:,1),centro(:,2),...
    [VV_ex(:,1); NaN; VV_in(:,1)], [VV_ex(:,2); NaN; VV_in(:,2)]);
ind_t_vess=find(ind_t_vess);
    triplot(two_d_elements(ind_t_vess,:),two_d_nodes(:,1),two_d_nodes(:,2),'y')
% Shell
plot(Shel_in(:,1),Shel_in(:,2),'k')
plot(Shel_ex(:,1),Shel_ex(:,2),'k')
[ind_t_shell,~]=find(inpolygon(centro(:,1),centro(:,2),Shel_ex(:,1),Shel_ex(:,2)));
[ind,~]=find(inpolygon(centro(:,1),centro(:,2),Shel_in(:,1), Shel_in(:,2)));
ind_t_shell=setdiff(ind_t_shell,ind);
    triplot(two_d_elements(ind_t_shell,:),two_d_nodes(:,1),two_d_nodes(:,2),'g')
% Saddle
load([input_Data ,'saddle.mat']);
ind_t_saddle=[];
for ii=1:4
    Points=saddle{ii};
    [Points] = fun_ordinapunti(Points);
    ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    ind_t_saddle=[ind_t_saddle; ind];
end
    triplot(two_d_elements(ind_t_saddle,:),two_d_nodes(:,1),two_d_nodes(:,2),'r')
% TSS
[ind_t_TSS,~]=find(inpolygon(centro(:,1),centro(:,2),TSS_ex(:,1),TSS_ex(:,2)));
[ind,~]=find(inpolygon(centro(:,1),centro(:,2),TSS_in(:,1), TSS_in(:,2)));
ind_t_TSS=setdiff(ind_t_TSS,ind);
ind_t_TSS=setdiff(ind_t_TSS,ind_t_saddle);
    triplot(two_d_elements(ind_t_TSS,:),two_d_nodes(:,1),two_d_nodes(:,2),'c')

% Conduttori (OH + FS)
% % ind_conduttori={};
% % for ii=1:40
% %     Points=cond{ii};
% %     ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
% %     ind_conduttori{ii}=ind;
% %     c=rand(3,1);
% %     triplot(two_d_elements(ind,:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
% % end
% % for ii=41:2:48
% %     Points=cond{ii};
% %     ind1=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
% %     Points=cond{ii+1};
% %     ind2=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
% %     ind_conduttori{ii}=[ind1; ind2];
% %     c=rand(3,1);
% %     triplot(two_d_elements([ind1; ind2],:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
% % end    
% % for ii=49:60
% %     Points=cond{ii};
% %     ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
% %     ind_conduttori{ii}=ind;
% %     c=rand(3,1);
% %     triplot(two_d_elements(ind,:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
% % end

ind_conduttori=[];
for ii=1:40
    Points=cond{ii};
    ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    ind_conduttori=[ind_conduttori; ind];
    c=rand(3,1);
    triplot(two_d_elements(ind,:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
end
for ii=41:2:48
    Points=cond{ii};
    ind1=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    Points=cond{ii+1};
    ind2=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    ind_conduttori=[ind_conduttori; ind1; ind2];
    c=rand(3,1);
    triplot(two_d_elements([ind1; ind2],:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
end    
for ii=49:60
    Points=cond{ii};
    ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    ind_conduttori=[ind_conduttori; ind];
    c=rand(3,1);
    triplot(two_d_elements(ind,:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
end

% Conduttori (Saddle)
load([input_Data ,'saddle.mat']);
ind_saddle={};
for ii=1:4
    Points=saddle{ii};
    [Points] = fun_ordinapunti(Points);
    ind=find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    ind_saddle{ii}=ind;
    c=rand(3,1);
    triplot(two_d_elements(ind,:),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)    
end

% nodi boundary vacuum vessel
ind_n_in=reshape(two_d_elements(ind_t_pla,:),numel(two_d_elements(ind_t_pla,:)),1);
ind_n_in=unique(ind_n_in);
ind_n=boundary(two_d_nodes(ind_n_in,1),two_d_nodes(ind_n_in,2),0.2);
ind_n=ind_n_in(ind_n);
ind_n_vess=(ind_n(1:end-1));
        
% % plot(two_d_nodes(ind_n_vess,1),two_d_nodes(ind_n_vess,2),'*c') 
        
%nodi nel boundary box
ind_bound=boundary(two_d_nodes(:,1),two_d_nodes(:,2),0.3);
ind_bound=ind_bound(1:end-1);
    plot(two_d_nodes(ind_bound,1),two_d_nodes(ind_bound,2),'k*');

%% Keyreg
load([input_Data ,'RFX_create_l_saddle_01.mat'],'keyreg');
load([input_Data ,'RFX_create_l_saddle_01_linear.mat']);
keyreg_create=keyreg;
meshData_create=meshData;

nt=size(two_d_elements,1);
keyreg=zeros(1,nt);
point.RR=meshData_create.n(:,1);
point.ZZ=meshData_create.n(:,2);
centro_coils=zeros(size(cond,2),2);
for ii=1:size(cond,2)
    centro_coils(ii,:)=sum(cond{ii})/size(cond{ii},1);
end
for j=1:14
    c=rand(3,1);
    flag_circuiti=find(KONNAX(j,:));
    ind_circuito=[];
    for k=1:numel(flag_circuiti)
        ind_t_subcirc=find(keyreg_create==flag_circuiti(k));
        ind_n_subcoil=unique(reshape(meshData_create.t(ind_t_subcirc,:),...
            numel(meshData_create.t(ind_t_subcirc,:)),1));
% %         triplot(meshData_create.t(ind_t_subcirc,:),point.RR,point.ZZ,'r')
% %         plot(point.RR(ind_n_subcoil),point.ZZ(ind_n_subcoil),'*','color',c)
        r_coil=meshData_create.n(ind_n_subcoil,1);
        z_coil=meshData_create.n(ind_n_subcoil,2);
        coils.rr=[min(r_coil); max(r_coil); max(r_coil); min(r_coil)];
        coils.zz=[min(z_coil); min(z_coil); max(z_coil); max(z_coil)];
        
% %         plot(coils.rr,coils.zz,'*k')
        %
        ind_tri=find(inpolygon(centro(:,1),centro(:,2),coils.rr,coils.zz));
        ind_tri=intersect(ind_tri,ind_conduttori);
% %         triplot(two_d_elements(ind_tri,1:3),two_d_nodes(:,1),two_d_nodes(:,2),'color',c)
        keyreg(ind_tri)=flag_circuiti(k);
% %         pause(0.0001)
    end
% %     pause
end

% Shell
keyreg(ind_t_shell)=117;
% Vessel
keyreg(ind_t_vess)=57;
% Saddle
for ii=1:4
    keyreg(ind_saddle{ii})=237+ii-1;
end
% TSS
keyreg(ind_t_TSS)=177;
% Plasma
keyreg(ind_t_pla)=-1;

%% Lati di bordo (tutto tranne il vuoto)
ind_vacuum=find(keyreg<=0);
ind_not_vacuum=setdiff(1:numel(keyreg),ind_vacuum);
triplot(two_d_elements(ind_not_vacuum,:),two_d_nodes(:,1),two_d_nodes(:,2),'k')
LATIBORDO=1;
meshData.t_choosen=two_d_elements(ind_not_vacuum,:);
[~,lati_bordo]=fun_lati(meshData,LATIBORDO);
figure
edges=pdemesh(two_d_nodes',lati_bordo',[]);
    set(edges,'color','k')
    set(edges,'linewidth',1)  
    
%% Midpoints (quadratic interpolation)
% % % Lati (tutto)
% % meshData.t_choosen=two_d_elements;
% % [lati,lati_bordo]=fun_lati(meshData,0);
% %     
% % midpoints=zeros(size(lati,1),2);
% % for ii=1:size(lati,1)
% %     midpoints(ii,:)=sum(two_d_nodes(lati(ii,:),:))*0.5;
% % end
% % two_d_nodes=[two_d_nodes; midpoints];
% % 
% % jj=size(two_d_nodes,1)+1;
% % midpoints=zeros(size(two_d_nodes));
% % midpoints_index=zeros(size(two_d_elements));
% % tic
% % for ii=1:size(two_d_elements,1)
% % % %     ii
% %     triangle=two_d_elements(ii,:);
% %     P=two_d_nodes(triangle,:);
% %     M=[0.5*(P(1,:)+P(2,:))
% %         0.5*(P(2,:)+P(3,:));...
% %         0.5*(P(3,:)+P(1,:))];
% %     a=ismember(M(1,:),midpoints,'rows');
% %     b=ismember(M(2,:),midpoints,'rows');
% %     c=ismember(M(3,:),midpoints,'rows');
% %     if ~a
% %         midpoints(ii,:)=M(1,:);
% %         mid(1)=jj;
% %         jj=jj+1;
% %     else
% %         mid(1)=size(two_d_nodes,1)+find(ismember(midpoints,M(1,:),'rows'));
% %     end
% %     if ~b
% %         midpoints(ii,:)=M(2,:);
% %         mid(2)=jj;
% %         jj=jj+1;
% %     else
% %         mid(2)=size(two_d_nodes,1)+find(ismember(midpoints,M(2,:),'rows'));
% %     end
% %     if ~c
% %         midpoints(ii,:)=M(3,:);
% %         mid(3)=jj;
% %         jj=jj+1;
% %     else
% %         mid(3)=size(two_d_nodes,1)+find(ismember(midpoints,M(3,:),'rows'));
% %     end
% %     midpoints_index(ii,:)=mid;
% % end
% % toc
% % two_d_elements=[two_d_elements midpoints_index];
% % two_d_nodes=[two_d_nodes; midpoints];
% % 
% % figure
% % pdemesh(two_d_nodes',[],[two_d_elements ones(size(two_d_elements,1),1)]');

%% Salvataggio dati (file meshData)
clear meshData
meshData.t=two_d_elements;
meshData.n=two_d_nodes;
meshData.e=lati_bordo;

%% Save Data
meshData.type=keyreg;
meshData.b=ind_bound;
meshData.vess=ind_n_vess;
meshData.type_explanation{1,1}='keyreg -1: plasma region';
meshData.type_explanation{2,1}='keyreg > 0: coils';
meshData.type_explanation{3,1}='keyreg = 117: shell';
meshData.type_explanation{4,1}='keyreg = 57: vessel';
meshData.type_explanation{5,1}='keyreg = 237: saddle coils';
meshData.type_explanation{6,1}='keyreg = 177: TSS';
save RFX_mesh_linear meshData

%% Midpoints dei lati 
% % %edges
% % tri_srt=two_d_elements;
% % ntri=size(two_d_elements,1);
% % lati=zeros(1,2);
% % lati_tot=zeros(3*size(two_d_elements,1),2);
% % for i=1:ntri
% %     triangle=tri_srt(i,:);
% %     lati_loc=[triangle(1) triangle(2);
% %         triangle(2) triangle(3);
% %         triangle(3) triangle(1)];
% %     if i==1
% %         lati=lati(2:end,:);
% %     end
% %     lati_tot(3*i-2:3*i,:)=lati_loc;
% % end
% % lati_tot=sort(lati_tot,2);
% % [lati,~,~]=unique(lati_tot,'rows');
% % midpoints=zeros(size(lati,1),2);
% % for ii=1:size(lati,1)
% %     ind_P1=lati(ii,1);
% %     ind_P2=lati(ii,2);
% %     P1=meshData.n(ind_P1,:);
% %     P2=meshData.n(ind_P2,:);
% %     midpoints(ii,:)=0.5*sum([P1; P2]);
% %     
% % end
% % % % plot(midpoints(:,1),midpoints(:,2),'r.') 
% % meshData.n=[meshData.n; midpoints];
% % 
% % for ii=1:ntri
% %     triangolo=two_d_elements(ii,:);
% %     P1=meshData.n(triangolo(1,:));
% %     P2=meshData.n(triangolo(2,:));
% %     P3=meshData.n(triangolo(3,:));
% %     M1=0.5*sum([P1; P2]);
% %     M2=0.5*sum([P2; P3]);
% %     M3=0.5*sum([P3; P1]);   
% % end
% % 
% % meshData.t=[meshData.t ind_lati_tri];
% % save RFX_mesh_quadratic meshData
% % 
% % figure
% % c=[0.5 0.5 0.5];
% % pdemesh(meshData.n',[],[meshData.t ones(size(meshData.t,1),1)]');
% % axis equal
% % axis([0.8 4.5 -1.5 1.5])
% % hold on

