function [meshData_loc] = fun_meshData_loc_4FRIDA(meshData,comp_domain)


%% Computational Domain (meshData_loc)

if nargin == 1
    comp_domain = meshData.n(meshData.ind_n_bc,:);
end

ind_t_comp=find(inpolygon(meshData.c_t(:,1),meshData.c_t(:,2),comp_domain(:,1),comp_domain(:,2)));
meshData_loc.ind_t_comp=ind_t_comp;

temp=meshData.t(ind_t_comp,:)';
ind_n_comp=unique(temp(:),'stable');
meshData_loc.ind_n_comp = ind_n_comp;

% % figure; hold on; axis equal
% % plot(meshData.n(ind_n_comp,1),meshData.n(ind_n_comp,2),'ro')


node_new=meshData.n(ind_n_comp,:);

if meshData.shape_functions_N_order == 1
    tri_new=zeros(length(ind_t_comp),3);
elseif meshData.shape_functions_N_order == 2
    tri_new=zeros(length(ind_t_comp),6);
end

tri_chosen=meshData.t(ind_t_comp,:);
for ii=1:length(ind_n_comp)
    [a,b]=find(tri_chosen==ind_n_comp(ii));
    for jj=1:length(a)
        tri_new(a(jj),b(jj))=ii;
    end
end

% % pdemesh(node_new',[],[tri_new(:,1:3)'; ones(1,size(tri_new,1))],'EdgeColor','c'); 

meshData_loc.t=tri_new;
meshData_loc.n=node_new;
meshData_loc.type=meshData.type(ind_t_comp);
meshData_loc.nn=size(meshData_loc.n,1);
meshData_loc.nt = size(meshData_loc.t,1);

meshData_loc.shape_functions_N_order = meshData.shape_functions_N_order;


P1=[meshData_loc.n(meshData_loc.t(:,1),1), ...
    meshData_loc.n(meshData_loc.t(:,1),2)];
P2=[meshData_loc.n(meshData_loc.t(:,2),1), ...
    meshData_loc.n(meshData_loc.t(:,2),2)];
P3=[meshData_loc.n(meshData_loc.t(:,3),1), ...
    meshData_loc.n(meshData_loc.t(:,3),2)];
meshData_loc.c_t=[sum(P1(:,1)+P2(:,1)+P3(:,1),2)/3, ...
    sum(P1(:,2)+P2(:,2)+P3(:,2),2)/3];

ind_t_loc_InFW=find(meshData_loc.type==-1)';
ind_n_loc_InFW=unique(meshData_loc.t(ind_t_loc_InFW,:),'stable');

meshData_loc.ind_t_InFW = ind_t_loc_InFW;
meshData_loc.ind_n_InFW = ind_n_loc_InFW;

% % meshData_loc.vess=boundary(meshData_loc.n(ind_n_loc_Invess,1),...
% %     meshData_loc.n(ind_n_loc_Invess,2),.00005);
% % meshData_loc.vess=ind_n_loc_Invess(meshData_loc.vess(1:end-1));

n_points_1 = size(comp_domain,1);

ind_n_loc_Onvess = zeros(n_points_1,1);
for ii = 1:n_points_1
    point_ii = comp_domain(ii,:);
    
    distance = sqrt((point_ii(1)-meshData_loc.n(:,1)).^2 + ...
       (point_ii(2)-meshData_loc.n(:,2)).^2);
    [~,ind_value] = min(distance);
    ind_n_loc_Onvess(ii) = ind_value;
           
end

% % meshData_loc.vess = ind_n_loc_Onvess;



% Area of primal and dual elements (local mesh)
[Area_tri_loc]=fun_AreaMesh_FEM(meshData_loc);
meshData_loc.Area_tri=Area_tri_loc;

% % ind_all2loc=zeros(meshData_loc.nn,1);
% % 
% % for ii=1:meshData_loc.nn
% %    P_loc_ii=meshData_loc.n(ii,:);
% %    dist_ii=sqrt((P_loc_ii(1)-meshData.n(:,1)).^2 + ...
% %        (P_loc_ii(2)-meshData.n(:,2)).^2);
% %    [aa,bb]=min(dist_ii);
% %    if aa == 0
% %        ind_all2loc(ii)=bb;
% %    else
% %        error('wrong point found')
% %    end 
% %    
% % end
% % 
% % meshData.ind_all2loc=ind_all2loc;
% % meshData_loc.ind_all2loc=ind_all2loc;

% % plot(meshData.n(ind_all2loc,1),meshData.n(ind_all2loc,2),'g*')

%% Plasma Domain and boundary conditions on local mesh

ind_FW = zeros(size(meshData.ind_n_FW,1),1);
for ii = 1:length(ind_FW)
    
    point_ii = meshData.n(meshData.ind_n_FW(ii),:);
    distance = sqrt((point_ii(1)-meshData_loc.n(:,1)).^2 + (point_ii(2)-meshData_loc.n(:,2)).^2);
    [~,bb] = min(distance);
    ind_FW(ii) = bb;
    
end

% % gamma = meshData_loc.n(ind_FW,:);
% % [~,ind_TH] = fun_ordinapunti(gamma);
% % ind_FW = ind_FW(ind_TH);

meshData_loc.ind_n_FW = ind_FW;


ind_n_bc = zeros(size(comp_domain,1),1);
for ii = 1:size(comp_domain,1)
    
    point_ii = comp_domain(ii,:);
    distance = sqrt((point_ii(1)-meshData_loc.n(:,1)).^2 + (point_ii(2)-meshData_loc.n(:,2)).^2);
    [~,bb] = min(distance);
    ind_n_bc(ii) = bb;
    
end

% % gamma = meshData_loc.n(ind_n_bc,:);
% % [~,ind_TH] = fun_ordinapunti(gamma);
% % ind_n_bc = ind_n_bc(ind_TH);

meshData_loc.ind_n_bc = ind_n_bc;


% % plot(meshData_loc.n(meshData_loc.vess,1),meshData_loc.n(meshData_loc.vess,2),'*m')
% % plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'or')
% % plot(meshData_loc.n(meshData_loc.ind_n_bc,1),meshData_loc.n(meshData_loc.ind_n_bc,2),'ob')


ind_D=setdiff(1:meshData_loc.nn,meshData_loc.ind_n_bc);
nn_D=length(ind_D);
nn_B=length(meshData_loc.ind_n_bc);

ind_B = meshData_loc.ind_n_bc;

if size(ind_D,1) == 1; ind_D = ind_D'; end
if size(ind_B,1) == 1; ind_D = ind_B'; end

meshData_loc.ind_D = ind_D;
meshData_loc.ind_B = ind_B;
meshData_loc.nn_D = nn_D;
meshData_loc.nn_B = nn_B;


% % meshData_loc.vers_n = meshData.vers_n;


%%
figure
subplot(1,2,1)
triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2),'k')
axis equal; hold on;
triplot(meshData.t(meshData.type == -1,1:3),meshData.n(:,1),meshData.n(:,2),'r')
title('meshData')
subplot(1,2,2)
triplot(meshData.t(:,1:3),meshData.n(:,1),meshData.n(:,2),'color', [.9 .9 .9])
axis equal; hold on;
triplot(meshData_loc.t(:,1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'k')
triplot(meshData_loc.t(meshData_loc.type == -1,1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'r')
plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'c*','LineWidth',2) 
plot(meshData_loc.n(meshData_loc.ind_n_bc,1),meshData_loc.n(meshData_loc.ind_n_bc,2),'og','LineWidth',2) 
title('meshData local')
% % quiver(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2), ...
% %     meshData_loc.vers_n(:,1),meshData_loc.vers_n(:,2))




































