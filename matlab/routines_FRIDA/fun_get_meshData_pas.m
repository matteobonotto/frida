function meshData_pas = fun_get_meshData_pas(meshData_ext)

%% Select only passive elements

ind_pas = meshData_ext.ind_pas;
if numel(ind_pas) == 1
    ind_t_pas =  find(meshData_ext.type == ind_pas);
else
    ind_t_pas = fun_vec_find(meshData_ext.type,ind_pas);
    % %     ind_t_pas =  find(sum(meshData_ext.type.' == (ind_pas)))';
end

ind_t_comp = ind_t_pas;

temp = meshData_ext.t(ind_t_comp,:)';
ind_n_comp=unique(temp(:),'stable');



node_new=meshData_ext.n(ind_n_comp,:);
N_order= meshData_ext.N_order;

if N_order == 1
    tri_new=zeros(length(ind_t_comp),3);
elseif N_order == 2
    tri_new=zeros(length(ind_t_comp),6);
end

tri_chosen = meshData_ext.t(ind_t_comp,:);
for ii=1:length(ind_n_comp)
    [a,b]=find(tri_chosen==ind_n_comp(ii));
    for jj=1:length(a)
        tri_new(a(jj),b(jj))=ii;
    end
end

if N_order == 1
    tri_new=tri_new(:,1:3);
end


% sort triangles
if size(ind_pas,1)>size(ind_pas,2)
    ind_t_map = repmat(ind_pas,1,2);
else
    ind_t_map = repmat(ind_pas.',1,2);
end

type_old_loc = meshData_ext.type(ind_t_comp);
type_new_loc = zeros(size(tri_new,1),1);

for ii = 1:size(ind_t_map,1)
    if ind_t_map(ii,1) == ind_t_map(ii,2)
        ind_sel = find((type_old_loc == (ind_t_map(ii,1):ind_t_map(ii,2))));
    else
        ind_sel = find(sum(type_old_loc.' == (ind_t_map(ii,1):ind_t_map(ii,2)).'))';
    end
    type_new_loc(ind_sel) = ii;
end

temp_node_new = zeros(size(node_new));
temp_tri_new  = zeros(size(tri_new));

kk = 1;
for ii = 1:size(ind_t_map,1)
    ind_nodes_ii = unique(tri_new(find(type_new_loc == ii),:));
    for jj = 1:numel(ind_nodes_ii)
        temp_node_new(kk,:) = node_new(ind_nodes_ii(jj),:);
        [ind_row,ind_col,~] = find(tri_new == ind_nodes_ii(jj));
        for hh = 1:numel(ind_row)
            temp_tri_new(ind_row(hh),ind_col(hh)) = kk;
        end
        kk = kk+1;
    end
end

tri_new = temp_tri_new;
node_new = temp_node_new;

meshData_pas.t = tri_new;
meshData_pas.n = node_new;
meshData_pas.type = meshData_ext.type(ind_t_comp);
meshData_pas.nn = size(meshData_pas.n,1);
meshData_pas.nt = size(meshData_pas.t,1);
meshData_pas.N_order = N_order;

