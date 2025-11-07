function [meshData_loc]=fun_UpdateNodes_fast(meshData_loc,solk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Update the nodes of the mesh that are inside the new boundary (plasma
%   nodes)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% compute gradient on nodes and triangle



%% Triangles and Nodes inside Boundary (geometrical)
Boundary=solk.Separatrix;
Nodes=meshData_loc.n;
ind_n_INpla=find(inpolygon(Nodes(:,1),Nodes(:,2),Boundary(:,1),Boundary(:,2)));

Baricentri = meshData_loc.c_t(meshData_loc.ind_t_InFW,:);
ind_t_INpla = inpolygon(Baricentri(:,1),Baricentri(:,2),Boundary(:,1),Boundary(:,2));
ind_t_INpla = meshData_loc.ind_t_InFW(ind_t_INpla);
% % ind_n_INpla_1=unique(meshData.t(ind_t_INpla,:));


%% output
meshData_loc.ind_t_INpla = ind_t_INpla;
meshData_loc.ind_n_INpla = ind_n_INpla;



%% Nodes on Boundary

% % nodes_pla=meshData_loc.n(ind_n_INpla,:);
% % ind_bound = boundary(nodes_pla(:,1),nodes_pla(:,2),.8);
% % ind_bound = unique(ind_bound);

% % meshData_loc.ind_n_ONpla = ind_n_INpla(ind_bound);



%% Nodes outside plasma
% % ind_n_NOTpla = setdiff(1:meshData.nn,[index.ind_n_INpla; index.ind_n_ONpla]);
ind_n_NOTpla = setdiff(1:meshData_loc.nn,ind_n_INpla);
meshData_loc.ind_n_NOTpla = ind_n_NOTpla;









end

