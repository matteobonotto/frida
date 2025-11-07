function [meshData_loc]=fun_UpdateNodes(meshData_loc,solk,SETTINGS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Update the nodes of the mesh that are inside the new boundary (plasma
%   nodes)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Triangles and Nodes inside Boundary (geometrical)
% % Boundary=solk.Separatrix;
% % Nodes=meshData_loc.n;
% % ind_n_INpla=find(inpolygon(Nodes(:,1),Nodes(:,2),Boundary(:,1),Boundary(:,2)));
% % 
% % ind_n_NOTpla = setdiff(1:meshData_loc.nn,ind_n_INpla);

Points     = meshData_loc.n;
psi_P      = solk.Psi;
Grad_psi_P = solk.Grad_nodes;
psi_B      = solk.Psi_B;
Point_Axis = [solk.Axis_RR solk.Axis_ZZ];
[ind_n_INpla,ind_n_NOTpla] = fun_find_points_INandOUT_separatrix(Points,psi_P,Grad_psi_P,psi_B,Point_Axis);

%%% Nodes on Boundary
if SETTINGS.RUN_FIXBOUNDARY == 1
    nodes_pla=meshData_loc.n(ind_n_INpla,:);
    ind_bound = boundary(nodes_pla(:,1),nodes_pla(:,2),.8);
    ind_bound = unique(ind_bound);
    ind_n_ONpla = ind_n_INpla(ind_bound);
else
    ind_n_ONpla = [];
end


% % Rplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,1))...
% %     max(meshData_loc.n(meshData_loc.ind_n_bc,1))]+[-0.3 0.3];
% % Zplot=[min(meshData_loc.n(meshData_loc.ind_n_bc,2))...
% %     max(meshData_loc.n(meshData_loc.ind_n_bc,2))]+[-0.3 0.3];
% % 
% % rr=meshData_loc.n(:,1);
% % zz=meshData_loc.n(:,2);
% % 
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % 
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % 
% % 
% %     PSI = griddata(rr,zz,solk0.Psi_MS,RR,ZZ);
% %     contour(RR,ZZ,PSI,50);



% % aa = tic;
% % for ii = 1:5000
% % [ind_n_in,ind_n_out] = fun_find_points_INandOUT_separatrix(Points,psi_P,Grad_psi_P,psi_B,Point_Axis);
% % end
% % toc(aa)/ii

%% Triangles and Nodes inside Boundary (geometrical)


% % Points     = meshData_loc.c_t;
% % psi_P      = solk.Psi_c_t;
% % Grad_psi_P = solk.Grad_c_t;
% % psi_B      = solk.Psi_B;
% % Point_Axis = [solk.Axis_RR solk.Axis_ZZ];
% % [ind_t_INpla,~] = fun_find_points_INandOUT_separatrix(Points,psi_P,Grad_psi_P,psi_B,Point_Axis);

qq = zeros(meshData_loc.nn,1);
qq(ind_n_INpla)= 1;
ind_t_INpla = find(sum(qq(meshData_loc.t),2)>0);

% % figure
% % triplot(meshData_loc.t(ind_t_INpla2,1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'k')
% % hold on; axis equal
% % triplot(meshData_loc.t(ind_t_INpla,1:3),meshData_loc.n(:,1),meshData_loc.n(:,2),'r')
% % plot(solk.Separatrix(:,1),solk.Separatrix(:,2),'k','linewidth',2)

% % Baricentri = meshData_loc.c_t(meshData_loc.ind_t_InFW,:);
% % ind_t_INpla2 = inpolygon(Baricentri(:,1),Baricentri(:,2),Boundary(:,1),Boundary(:,2));
% % ind_t_INpla2 = meshData_loc.ind_t_InFW(ind_t_INpla2);

% % figure
% % plot(meshData_loc.c_t(:,1),meshData_loc.c_t(:,2),'.k')
% % hold on; axis equal;
% % plot(meshData_loc.c_t(ind_t_INpla,1),meshData_loc.c_t(ind_t_INpla,2),'or')
% % plot(solk.Separatrix(:,1),solk.Separatrix(:,2),'k','linewidth',2)
% % plot(meshData_loc.n(ind_n_INpla,1),meshData_loc.n(ind_n_INpla,2),'dg')
% % 
% % figure
% % plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),psi_P,'o')


%% Gauss points inside plasma (psi>psi_b & grad(psi)*vec_Or<0)

Points     = meshData_loc.nodes_Gauss_pla;
psi_P      = solk.Psi_Gauss;
Grad_psi_P = solk.Grad_Gauss;
psi_B      = solk.Psi_B;
Point_Axis = [solk.Axis_RR solk.Axis_ZZ];
[ind_G_in,ind_G_out] = fun_find_points_INandOUT_separatrix(Points,psi_P,Grad_psi_P,psi_B,Point_Axis);

% % aa = tic;
% % for ii = 1:1000
% % [ind_G_in,ind_G_out] = fun_find_points_INandOUT_separatrix(Points,psi_P,Grad_psi_P,psi_B,Point_Axis);
% % end
% % toc(aa)/ii
% % 
% % ind_sel_1 = find(solk.Psi_Gauss >= solk.Psi_B);
% % % % Psi_Gauss_sel = solk.Psi_Gauss(ind_sel_1);
% % 
% % vec_Or = meshData_loc.nodes_Gauss_pla - [solk.Axis_RR solk.Axis_ZZ];
% % 
% % prod_scal = sum(solk.Grad_Gauss.*vec_Or,2);
% % ind_sel_2 = find(prod_scal <= 0);
% % 
% % ind_sel_final = intersect(ind_sel_1,ind_sel_2);
% % 
% % figure
% % plot(meshData_loc.nodes_Gauss_pla(ind_sel,1),meshData_loc.nodes_Gauss_pla(ind_sel,2),'o')
% % plot(meshData_loc.nodes_Gauss_pla(ind_sel_2,1),meshData_loc.nodes_Gauss_pla(ind_sel_2,2),'o')
% % plot(meshData_loc.nodes_Gauss_pla(ind_sel_final,1),meshData_loc.nodes_Gauss_pla(ind_sel_final,2),'o')
% % 
% % ind_sel_1 = (solk.Psi_Gauss >= solk.Psi_B);
% % % % Psi_Gauss_sel = solk.Psi_Gauss(ind_sel_1);
% % 
% % vec_Or = meshData_loc.nodes_Gauss_pla - [solk.Axis_RR solk.Axis_ZZ];
% % 
% % prod_scal = sum(solk.Grad_Gauss.*vec_Or,2);
% % ind_sel_2 = (prod_scal <= 0);
% % 
% % ind_G_in = (ind_sel_1 & ind_sel_2);
% % ind_G_out = (true(meshData_loc.n_Gauss*meshData_loc.nt,1) & ~ind_G_in);
% % 
% % figure
% % plot(meshData_loc.nodes_Gauss_pla(:,1),meshData_loc.nodes_Gauss_pla(:,2),'*k')
% % hold on; axis equal;
% % plot(meshData_loc.nodes_Gauss_pla(ind_G_in,1),meshData_loc.nodes_Gauss_pla(ind_G_in,2),'or')
% % plot(meshData_loc.nodes_Gauss_pla(ind_G_out,1),meshData_loc.nodes_Gauss_pla(ind_G_out,2),'dg')
% % plot(solk.Separatrix(:,1),solk.Separatrix(:,2),'k','linewidth',2)

%% output
meshData_loc.ind_t_INpla  = ind_t_INpla;
meshData_loc.ind_n_INpla  = ind_n_INpla;
meshData_loc.ind_n_ONpla  = ind_n_ONpla;
meshData_loc.ind_n_NOTpla = ind_n_NOTpla;
meshData_loc.ind_G_in     = ind_G_in;
meshData_loc.ind_G_out    = ind_G_out;








end

