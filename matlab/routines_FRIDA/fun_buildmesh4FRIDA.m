function [meshData] = fun_buildmesh4FRIDA(...
    FW, ...
    n_FW, ...
    bc, ...
    n_bc, ...
    factor_bc, ...
    Conductors_geo, ...
    Conductors_index, ...
    spacing, ...
    N_order, ...
    dir_FRIDA, ...
    MEX_OPT)

%% Domains for meshing

figure
plot(FW(:,1),FW(:,2),'or-')
axis equal; hold on;

% resampling bc (if n_FW is defined)
if ~isempty(n_FW)
    [FW_r,FW_z]=equispaced_stable(FW(:,1),FW(:,2),n_FW);
    % %     plot(FW_r,FW_z,'*-k')
    FW = [FW_r(1:end) FW_z(1:end)];
end


if isempty(bc)
    
    % boundary of FW
    ind = boundary(FW(:,1),FW(:,2),0);
    bc = FW(ind,:);
    plot(bc(:,1),bc(:,2),'o-')
% %     [bc,~] = fun_ordinapunti(bc);
    
    % % figure
    % % plot(FW(:,1),FW(:,2),'o-')
    % % axis equal; hold on
    % % plot(bc(:,1),bc(:,2),'o-')
    
    bc = unique(bc,'rows');
        
    [bc,~] = fun_ordinapunti(bc);
    [RR_EQUI,ZZ_EQUI]=equispaced_stable(bc(:,1),bc(:,2),n_bc);
    bc = [RR_EQUI ZZ_EQUI];
    
    U_t = zeros(size(bc));
    U_n = zeros(size(bc));
    nfil=length(bc(:,1));
    rstart=bc(:,1);
    zstart=bc(:,2);
    for kpoint=1:nfil
        kprev=kpoint-1; 
        if(kprev<1) kprev=nfil; end
        
        kfol=kpoint+1; 
        if(kfol>nfil) kfol=1; end
        
        t=[rstart(kfol)-rstart(kprev); zstart(kfol)-zstart(kprev)];
        t=t/norm(t);
        U_t(kpoint,:)=t';
        
        U_n(kpoint,:)=[t(2) -t(1)]; 
    end
        
    bc = bc + factor_bc*U_n;
% %     plot(bc(:,1),bc(:,2),'o-')

else

    [bc_r,bc_z]=equispaced_stable(bc(:,1),bc(:,2),n_bc);
    bc = [bc_r(1:end) bc_z(1:end)];

end


%% 2D meshing with GMSH

GEOMETRY.shape{1} = FW;
GEOMETRY.shape{2} = bc;
for ii=1:numel(Conductors_geo)
    GEOMETRY.shape{ii+2}=Conductors_geo{ii};
end

GEOMETRY.spacing{1} = spacing(1);
GEOMETRY.spacing{2} = spacing(2);
for ii=1:numel(Conductors_geo)
    GEOMETRY.spacing{ii+2}=spacing(ii+2);
end

GEOMETRY.surface{1} = 1;
GEOMETRY.surface{2} = [2 1];
for ii=1:numel(Conductors_geo)
    GEOMETRY.surface{ii+2}=ii+2;
end

mesh_name = 'temp_mesh';
[p,e,t]=autoMESH_complex(mesh_name,GEOMETRY,N_order,dir_FRIDA);


figure
triplot(t(:,1:3),p(:,1),p(:,2),'k')
axis equal; hold on;
% % plot(FW(:,1),FW(:,2),'o')

%%
meshData.n = p;
meshData.e = e;
meshData.t = t;

nn = size(meshData.n,1);
nt = size(meshData.t,1);
ne = size(meshData.e,1);

meshData.nn = nn;
meshData.nt = nt;
meshData.ne = ne;

%% Find index insidee various domains
type=zeros(meshData.nt,1);

%baricentro triangoli
if N_order == 1
    P1=[meshData.n(meshData.t(:,1),1), ...
        meshData.n(meshData.t(:,1),2)];
    P2=[meshData.n(meshData.t(:,2),1), ...
        meshData.n(meshData.t(:,2),2)];
    P3=[meshData.n(meshData.t(:,3),1), ...
        meshData.n(meshData.t(:,3),2)];
    
    centro=[sum(P1(:,1)+P2(:,1)+P3(:,1),2)/3, ...
        sum(P1(:,2)+P2(:,2)+P3(:,2),2)/3];
    
elseif N_order == 2
    P1=[meshData.n(meshData.t(:,1),1), ...
        meshData.n(meshData.t(:,1),2)];
    P2=[meshData.n(meshData.t(:,2),1), ...
        meshData.n(meshData.t(:,2),2)];
    P3=[meshData.n(meshData.t(:,3),1), ...
        meshData.n(meshData.t(:,3),2)];
    P4=[meshData.n(meshData.t(:,4),1), ...
        meshData.n(meshData.t(:,4),2)];
    P5=[meshData.n(meshData.t(:,5),1), ...
        meshData.n(meshData.t(:,5),2)];
    P6=[meshData.n(meshData.t(:,6),1), ...
        meshData.n(meshData.t(:,6),2)];
    
    centro=[sum(P1(:,1)+P2(:,1)+P3(:,1)+P4(:,1)+P5(:,1)+P6(:,1),2)/6, ...
        sum(P1(:,2)+P2(:,2)+P3(:,2)+P4(:,2)+P5(:,2)+P6(:,2),2)/6];

end


% plasma
[ind_t_pla,~]=inpolygon(centro(:,1),centro(:,2),FW(:,1),FW(:,2));
ind_t_pla=find(ind_t_pla);
type(ind_t_pla)=-1;

% Conductors
for ii=1:size(Conductors_geo,2)
    
    Points = Conductors_geo{ii};
% %     ind = find(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2)));
    type(inpolygon(centro(:,1),centro(:,2),Points(:,1), Points(:,2))) = Conductors_index(ii);
% % %     Conductors_index(ii)
% % %     triplot(meshData.t(type == Conductors_index(ii),1:3),meshData.n(:,1),meshData.n(:,2),'c'); 
% % %     hold on; axis equal
% % %     pause
end

triplot(t(ind_t_pla,1:3),p(:,1),p(:,2),'r')
% % triplot(t(type == 57,1:3),p(:,1),p(:,2),'r')

%% index of FW nodes and bc nodes
% FW nodes
ind_n_pla = meshData.t(type == -1,:);
ind_n_pla = ind_n_pla(:);

ind_n_vac = meshData.t(type == 0,:);
ind_n_vac = ind_n_vac(:);

ind_n_FW = intersect(ind_n_pla,ind_n_vac);


% bc nodes
ind_t_vac = find(type == 0);
ind_n_in=reshape(meshData.t(ind_t_vac,:),numel(meshData.t(ind_t_vac,:)),1);
ind_n_in=unique(ind_n_in);
ind_n=boundary(meshData.n(ind_n_in,1),meshData.n(ind_n_in,2),0.2);
ind_n=ind_n_in(ind_n);
ind_n_bc=(ind_n(1:end-1));


%
plot(meshData.n(ind_n_FW,1),meshData.n(ind_n_FW,2),'c*','LineWidth',2) 
plot(meshData.n(ind_n_bc,1),meshData.n(ind_n_bc,2),'og','LineWidth',2) 



%% vercots normal to FW (for limiter evaluation)


[FW_r,FW_z]=equispaced_stable(FW(:,1),FW(:,2),10000);

FW_new = [FW_r FW_z];

FW_mid = .5*(FW_new + FW_new([2:end 1],:));

vers_t = FW_mid - FW_mid([2:end 1],:);
vers_n = [-vers_t(:,2) vers_t(:,1)];

vers_n_r = scatteredInterpolant(FW_new(:,1),FW_new(:,2),vers_n(:,1));
vers_n_z = scatteredInterpolant(FW_new(:,1),FW_new(:,2),vers_n(:,2));

vers_n = [vers_n_r(meshData.n(ind_n_FW,:)) vers_n_z(meshData.n(ind_n_FW,:))];
vers_n = vers_n./fun_vecnorm(vers_n);

quiver(meshData.n(ind_n_FW,1),meshData.n(ind_n_FW,2),vers_n(:,1),vers_n(:,2))


meshData.vers_n = vers_n;

% % plot(FW_mid(:,1),FW_mid(:,2),'o')
% % quiver(FW_new(:,1),FW_new(:,2),vers_t(:,1),vers_t(:,2))
% % quiver(FW_new(:,1),FW_new(:,2),vers_n(:,1),vers_n(:,2))

% % RR = FW_new(:,1);
% % ZZ = FW_new(:,2);
% % 
% % DIST=sqrt((RR(2:end)-RR(1:end-1)).^2+(ZZ(2:end)-ZZ(1:end-1)).^2);
% % LENGTH = sum(DIST);
% % 
% % s_coo_ref=zeros(length(RR),1);
% % for ii=2:length(RR)
% %     s_coo_ref(ii)=s_coo_ref(ii-1)+DIST(ii-1);
% % end
% % 
% % 
% % RR = meshData.n(ind_n_FW,1);
% % ZZ = meshData.n(ind_n_FW,1);
% % DIST=sqrt((RR(2:end)-RR(1:end-1)).^2+(ZZ(2:end)-ZZ(1:end-1)).^2);
% % 
% % s_coo_FW=zeros(length(RR),1);
% % for ii=2:length(RR)
% %     s_coo_ref(ii)=s_coo_ref(ii-1)+DIST(ii-1);
% % end



%% interface edges

LATIBORDO=1;
t_choosen=meshData.t(type>=0,:);

if MEX_OPT 
    [~,e_interface]=fun_lati_fast_mex(t_choosen,LATIBORDO);
else
    [~,e_interface]=fun_lati_fast(t_choosen,LATIBORDO);
end

meshData.e_interface=e_interface;

% % figure
% % edges=pdemesh(meshData.n',e_interface',[]);
% % set(edges,'color','k')
% % set(edges,'linewidth',1)


%% Output (meshData)
meshData.type = type;
meshData.ind_n_FW=ind_n_FW;
meshData.ind_n_bc=ind_n_bc;
meshData.shape_functions_N_order=N_order;
meshData.c_t = centro;
meshData.FW= meshData.n(meshData.ind_n_FW,:);
meshData.bc= meshData.n(meshData.ind_n_bc,:);
meshData.type_explanation{1,1}='type == 0: vacuum';
meshData.type_explanation{2,1}='type == -1: plasma';
meshData.type_explanation{3,1}='type > 0: coils';

% % save([mesh_name '.mat'], 'meshData');































