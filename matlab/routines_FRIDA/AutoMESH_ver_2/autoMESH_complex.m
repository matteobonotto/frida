function [p,e,t]=autoMESH_complex(filename,GEOMETRY,order,dir_FRIDA)

%% WRITEB INPUT GMSH .geo FILE
writeGMSH(filename,GEOMETRY);

%% MESH THE GEOMETRY USING GMSH;

if ismac
    command_dos = sprintf('gmsh -2 %s.geo', filename);
else
    command_dos = sprintf('%s/AutoMESH_ver_2/gmsh -2 %s.geo', dir_FRIDA, filename);
end

dos(command_dos);
% save('mesh_geo', 'p','t')


%% READ MESH
[p,e,t]=readGMSH(filename);

if order == 2
    
    [p, e, t] = convertMeshToSecondOrder(p', e', t');
    
    p = p';
    e = e';
    t = t';
    
    t = t(:,1:6);
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function [ fileID ] = writeGMSH(filename,GEOMETRY)


shape = GEOMETRY.shape;
spacing = GEOMETRY.spacing;
surface = GEOMETRY.surface;

for ii=1:size(shape,2)
    refinement{ii}=spacing{ii}*ones(size(shape{ii},1),1);
end

%%

n_point = 1;
n_line  = 1;
n_loop  = 1;
n_surf  = 1;

fileID = fopen([filename '.geo'],'w');

for ii=1:size(shape,2)
    element=shape{ii};
    spacing=refinement{ii};
    for jj = 1:size(element,1)
        fprintf(fileID,'%s  %d %s %8f %s %8f %s %8f %s \n',...
            char('Point('),...
            n_point,...
            char(')={'),...
            element(jj,1),...
            ',',...
            element(jj,2),...
            ',0,',...
            spacing(jj),...
            '};');
        n_point=n_point+1;
    end
    fprintf(fileID,'%6s \n',' ');
    
    % LINES
    for jj = 1:size(element,1)
        if jj<size(element,1)
            fprintf(fileID,'%s  %d %s %d %s %d %s \n',...
                char('Line('),...
                n_line,...
                char(')={'),...
                n_line,...
                ',',...
                n_line+1,...
                '};');
        elseif jj==size(element,1)
            fprintf(fileID,'%s  %d %s %d %s %d %s \n',...
                char('Line('),...
                n_line,...
                char(')={'),...
                n_line,...
                ',',...
                n_line-size(element,1)+1,...
                '};');
        end
        n_line=n_line+1;
    end
    
    fprintf(fileID,'%1s \n',' ');
    
    % LINE LOOP
    fprintf(fileID,'Line Loop(%i)={',n_loop);
    for jj = 1:size(element,1)
        if jj<size(element,1)
            fprintf(fileID,'%d %s\n',...
                n_line-size(element,1)+jj-1,...
                ',');
        elseif jj==size(element,1)
            fprintf(fileID,'%d %s\n',...
                n_line-size(element,1)+jj-1, ...
                '};');
        end
    end
    n_loop = n_loop+1;

    fprintf(fileID,'%1s \n',' ');
    
    % PLANE SURFACE
    tmp_surf = surface{ii};
    
    tmp_string = 'Plane Surface(%i)={';
    for jj = 1:numel(tmp_surf)
        tmp_string = [tmp_string '%i,'];
    end
    tmp_string = tmp_string(1:end-1);
    tmp_string = [tmp_string '};'];
    
    fprintf(fileID,tmp_string,n_surf, tmp_surf);
    n_surf = n_surf+1;
    
    fprintf(fileID,'%1s \n\n',' ');

    
end

fclose(fileID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,e,t]=readGMSH(filename)

file=([filename '.msh']);

% no of nodes is mentioned in 5th row and first column
N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

node_id     = dlmread(file,'',[5 0 4+N_n 0]);
nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);

%------- 2D Geometry

p = nodes(:,1:2);
elem_type   = elements(:,2);

%% NODES
two_ind = 1;
for i = 1:N_e
    if(elem_type(i) ~= 2)
        two_ind = two_ind+1;
    end
end

%% ELEMENTS
t(1:N_e-two_ind,1:3) = 0;
k = 1;
for i = two_ind:N_e
    t(k,1:3) = elements(i,6:8);
    k = k+1;
end

%% EDGES
ntri=size(t,1);
lati=zeros(1,2);
lati_tot=zeros(3*size(t,1),2);
for i=1:ntri
    triangle=t(i,:);
    lati_loc=[triangle(1) triangle(2);
        triangle(2) triangle(3);
        triangle(3) triangle(1)];
    if i==1
        lati=lati(2:end,:);
    end
    lati_tot(3*i-2:3*i,:)=lati_loc;
end
lati_tot=sort(lati_tot,2);
[e,~,~]=unique(lati_tot,'rows');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p, e, t, nVnodes, nPnodes, indices, fluidIndices] = convertMeshToSecondOrder(p, e, t, varargin)

% convertMeshToSecondOrder - Convert linear grid to second order grid.
%
% This QuickerSim CFD Toolbox function converts a first order linear mesh
% read from an external file (e.g. Gmsh etc.) to the second order mesh
% which is needed for fluid flow simulation with P2-P1 elements (so called
% Taylor-Hood pair) with second order shape functions for velocity
% discretization and linear shape functions for pressure approximation.
%
% [p, e, t, nVnodes, nPnodes, indices] = convertMeshToSecondOrder(p, e, t)
% [p, e, t, nVnodes, nPnodes, indices] =
%                   convertMeshToSecondOrder(p, e, t, fluidDomainIds);
%
% Input arguments:
% p - array of nodal coordinates in the linear mesh.
% e - array of boundary elements in the linear mesh.
% t - array of linear domain elements (triangles/tetrahedra).
% fluidDomainIds - additional argument with ids of subdomains for coupled
%                  heat transfer between solid and fluid supported in the
%                  full toolbox version(for examples see e.g.
%                  genericHeatFluidSolidSolver2D).
%
% Output arguments:
% p, e, t - arrays defining mesh - see help for importMeshGmsh function for
%     detailed description of mesh PET format.
% p - as above with additional second order nodes in the middle of each
%     edge of every element.
% e - as above with added indices of new nodes. For a linear 2-D mesh this
%     array has 7 rows, for a second order grid the number of rows equals
%     8 - the last row stores the indices of the new nodes added in the
%     middle of each edge. For the 3-D case number of rows in e is either 6
%     or 9 for linear and second order grid, respectively. For details of
%     node numbering refer to documentation.
% t - as above. For a 2-D linear mesh this array stores indices of element
%     nodes in rows 1 to 3 and the id of the domain in the last 4-th row. In
%     case of second order grid the corner nodes of each triangle are still
%     stored in rows 1 to 3 and rows 4 to 6 store correspondingly indices of
%     new second order nodes in the following order: node 4 (placed between
%     original 1st and 2nd node), node 5 (placed between original 2nd and 3rd
%     node) and node 6 (placed between original 3rd and 1st node). In 3-D
%     number of rows in t equals 5 for a linear mesh (4 element nodes and
%     domain id) or 11 for second order mesh (10 ids of element nodes and
%     domain id in the last row). For details of node numbering refer to
%     documentation.
% nVnodes - total number of nodes associated with velocity unknowns (so called
%     velocity nodes - i.e. all nodes in the second order mesh).
% nPnodes - total number of nodes associated with pressure unknowns (so
%     called pressure nodes - i.e. all nodes of the original linear grid).
% indices - a structure containing for each node indices of the x-velocity,
%     y-velocity, z-velocity (only for 3-D case) and pressure unknown
%     in the global solution vector u.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also IMPORTMESHGMSH, GENERATEINDICES2D.

dim = size(p,1);

nnodes = size(p,2);
nelements = size(t,2);

Edges = cell(nnodes,1);
AddNodes = cell(nnodes,1);

for i = 1:nnodes
    Edges{i} = i;
    AddNodes{i} = i;
end

if(dim == 2)
    t = [t; zeros(3, nelements)];
    t(7,:) = t(4,:);
else % dim == 3
    t = [t; zeros(6, nelements)];
    t(11,:) = t(5,:);
end
el = nnodes+1;

if(dim == 2)
    for i = 1:nelements
        for j = 1:3
            n1 = t(j,i);
            n2 = t(mod(j,3)+1,i);
            nloc = sort([n1 n2]);
            n1 = nloc(1);
            n2 = nloc(2);
            
            ind = find(Edges{n1}==n2);
            
            if(ind)
                t(3+j,i) = AddNodes{n1}(ind);
            else
                Edges{n1}(end+1) = n2;
                AddNodes{n1}(end+1) = el;
                t(3+j,i) = el;
                el = el+1;
            end
        end
    end
else % dim == 3
    edgeNodes = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
    for i = 1:nelements
        for j = 1:6
            n1 = t(edgeNodes(j,1),i);
            n2 = t(edgeNodes(j,2),i);
            nloc = sort([n1 n2]);
            n1 = nloc(1);
            n2 = nloc(2);
            
            ind = find(Edges{n1}==n2);
            
            if(ind)
                t(4+j,i) = AddNodes{n1}(ind);
            else
                Edges{n1}(end+1) = n2;
                AddNodes{n1}(end+1) = el;
                t(4+j,i) = el;
                el = el+1;
            end
        end
    end
end

el = el-1;

przyrost = el-nnodes;

if(dim == 2)
    p = [p zeros(2,przyrost)];
    
    for i = 1:nelements
        for j = 1:3
            n1 = t(j,i);
            n2 = t(mod(j,3)+1,i);
            
            p(:,t(j+3,i)) = 0.5*(p(:,n1)+p(:,n2));
        end
    end
else
    p = [p zeros(3,przyrost)];
    
    edgeNodes = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
    for i = 1:nelements
        for j = 1:6
            n1 = t(edgeNodes(j,1),i);
            n2 = t(edgeNodes(j,2),i);
            
            p(:,t(j+4,i)) = 0.5*(p(:,n1)+p(:,n2));
        end
    end
end


% Enrich edges
if(dim == 2)
    nEdges = size(e,2);
    e = [e; zeros(1,nEdges)];
    
    for i = 1:nEdges
        n1 = e(1,i);
        n2 = e(2,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        ind = Edges{n1} == n2;
        e(8,i) = AddNodes{n1}(ind);
    end
    
    nVnodes = size(p,2);
    
    %     Renumber nodes to reduce problem matrix bandwidth
    % Assemble matrix exemplifying sparsity pattern
    %     kI = reshape(repmat(reshape(t(1:6, :), 1, 6*size(t, 2)), 6, 1), 36*size(t, 2), 1);
    %     kJ = reshape(repmat(t(1:6, :), 6, 1), 36*size(t, 2), 1);
    %     kV = ones(36*size(t, 2), 1);
    %     rrd = 1:size(p, 2);
    %     rrd(symrcm(sparse(kI, kJ, kV))) = rrd;
    %     temp = reshape(t(1:6, :), 6*size(t, 2), 1);
    %     t(1:6, :) = reshape(rrd(temp), 6, size(t, 2));
    %     temp = reshape(e([1:2, 8], :), 3*size(e, 2), 1);
    %     e([1:2, 8], :) = reshape(rrd(temp), 3, size(e, 2));
    %     p(:, rrd) = p;
    %
    nPnodes = max(max(t(1:3,:)));
    %     kI = reshape(repmat(reshape(t(4:6, :), 1, 3*size(t, 2)), 3, 1), 9*size(t, 2), 1) - nPnodes;
    %     kJ = reshape(repmat(t(4:6, :), 3, 1), 9*size(t, 2), 1) - nPnodes;
    %     kV = ones(9*size(t, 2), 1);
    %     rrd = 1:size(p, 2) - nPnodes;
    %     rrd(symrcm(sparse(kI, kJ, kV))) = rrd;
    %     temp = reshape(t(4:6, :), 3*size(t, 2), 1);
    %     t(4:6, :) = reshape(rrd(temp - nPnodes) + nPnodes, 3, size(t, 2));
    %     temp = reshape(e(8, :), size(e, 2), 1);
    %     e(8, :) = reshape(rrd(temp - nPnodes) + nPnodes, 1, size(e, 2));
    %     p(:, rrd + nPnodes) = p(:, nPnodes + 1:size(p, 2));
    %
    
    indices = generateIndices2D(p,t);
    
    if(~isempty(varargin))
        fluidDomainIds = varargin{1};
        fluidElements = ismember(t(7, :), fluidDomainIds);
        fluidNodes = unique(t(1:6,fluidElements));
        fluidNodesP = fluidNodes(fluidNodes <= nPnodes);
        fluidIndices = false(size(indices.indu));
        fluidIndices([fluidNodes(:); fluidNodes(:) + nVnodes; fluidNodesP(:) + 2*nVnodes]) = true;
        %         fluidIndices = [indices.indu(fluidNodes) indices.indv(fluidNodes) indices.indp(fluidNodesP)];
    end
else % dim == 3
    nEdges = size(e,2);
    e = [e(1:3,:); zeros(3,nEdges); e(4:6,:)];
    
    for i = 1:nEdges
        % Popraw cialo tej petli
        n1 = e(1,i);
        n2 = e(2,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        e(4,i) = AddNodes{n1}(Edges{n1}==n2);
        
        n1 = e(2,i);
        n2 = e(3,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        e(5,i) = AddNodes{n1}(Edges{n1}==n2);
        
        n1 = e(1,i);
        n2 = e(3,i);
        nloc = sort([n1 n2]);
        n1 = nloc(1);
        n2 = nloc(2);
        e(6,i) = AddNodes{n1}(Edges{n1}==n2);
    end
    
    nPnodes = max(max(t(1:4,:)));
    %     kI = reshape(repmat(reshape(t(5:10, :), 1, 6*size(t, 2)), 6, 1), 36*size(t, 2), 1) - nPnodes;
    %     kJ = reshape(repmat(t(5:10, :), 6, 1), 36*size(t, 2), 1) - nPnodes;
    %     kV = ones(36*size(t, 2), 1);
    %     rrd = 1:size(p, 2) - nPnodes;
    %     rrd(symrcm(sparse(kI, kJ, kV))) = rrd;
    %     temp = reshape(t(5:10, :), 6*size(t, 2), 1);
    %     t(5:10, :) = reshape(rrd(temp - nPnodes) + nPnodes, 6, size(t, 2));
    %     temp = reshape(e(4:6, :), 3*size(e, 2), 1);
    %     e(4:6, :) = reshape(rrd(temp - nPnodes) + nPnodes, 3, size(e, 2));
    %     p(:, rrd + nPnodes) = p(:, nPnodes + 1:size(p, 2));
    nVnodes = size(p,2);
    
    indices = generateIndices(p, t);
    
    if(~isempty(varargin))
        fluidDomainIds = varargin{1};
        fluidElements = ismember(t(11, :), fluidDomainIds);
        fluidNodes = unique(t(1:10, fluidElements));
        fluidNodesP = fluidNodes(fluidNodes <= nPnodes);
        fluidIndices = false(size(indices.indu));
        fluidIndices([fluidNodes(:); fluidNodes(:) + nVnodes; fluidNodes(:) + 2*nVnodes; fluidNodesP(:) + 3*nVnodes]) = true;
        %         fluidIndices = [indices.indu(fluidNodes) indices.indv(fluidNodes) indices.indw(fluidNodes) indices.indp(fluidNodesP)];
    end
    
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

function [indices] = generateIndices2D(p, t)

% generateIndices2D - Generate global ids of pressure and velocity unknowns
%
% This QuickerSim CFD Toolbox function generates indices for easy
% manipulation of the solution vector u.
%
% indices = generateIndices2D(p,t);
%
% Input arguments:
% p - array of nodal coordinates (for detailed description see help of the
%     importMeshGmsh function).
% t - array of finite elements (for detailed description see help for the
%     importMeshGmsh function).
%
% Output arguments:
% indices - a Matlab structure with fiels:
%           indices.indu - with indices of x-velocity unknowns,
%           indices.indv - with indices of y-velocity unknowns,
%           indices.indp - with indices of pressure unknowns.
% The general structure of the solution vector u in this Toolbox is that
% the first 1 to nVnodes (where nVnodes stands for total number of velocity
% nodes - see help for convertMeshToSecondOrder for detailed description of
% nVnodes) elements in the solution vector u correspond to x-velocity
% values in the subsequent nodes of the mesh, the next elements in u, this
% means u(nVnodes+1:2*nVnodes) contain values of the y-velocity for each
% mesh node and then the last elements, this means u(2*nVnodes+1:end)
% contain pressure values for each mesh node (note here that only first
% order nodes contain pressure solution - do not refer to mid-edge nodes
% when requesting pressure values, since this will cause an error with
% false indices in indices.indp vector). Function generatePressureData
% might be helpful.
%
% Examples:
% 1. Extract y-velocity value for the 179th node in the mesh:
%       indices = generateIndices2D(p,t);
%       v179 = u(indices.indv(179));
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: CONVERTMESHTOSECONDORDER, EXPORTTOGMSH2D, GENERATEPRESSUREDATA,
%           IMPORTMESHGMSH.

% Pressure nodes
pNodesNumber = max(max(t(1:3, :)));

% Velocity nodes
vNodesNumber = size(p, 2);

% Save indices as logical vectors
indu = false(vNodesNumber*2 + pNodesNumber, 1);
indv = indu;
indp = indu;

indu(1:vNodesNumber) = true;
indv(vNodesNumber + 1:2*vNodesNumber) = true;
indp((1:pNodesNumber) + 2*vNodesNumber) = true;

indices = struct('indp', indp, 'indu', indu, 'indv', indv);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

function [indices] = generateIndices(p, t)

% generateIndices - Generate global ids of pressure and velocity unknowns
%
% This QuickerSim CFD Toolbox function generates indices for easy
% manipulation of the solution vector u.
%
% indices = generateIndices(p,t);
%
% Input arguments:
% p - array of nodal coordinates (for detailed description see help of the
%     importMeshGmsh function).
% t - array of finite elements (for detailed description see help for the
%     importMeshGmsh function).
%
% Output arguments:
% indices - a Matlab structure with fiels:
%           indices.indu - with indices of x-velocity unknowns,
%           indices.indv - with indices of y-velocity unknowns,
%           indices.indw - with indices of z-velocity unknowns,
%           indices.indp - with indices of pressure unknowns.
% The general structure of the solution vector u in this Toolbox is that
% the first 1 to nVnodes (where nVnodes stands for total number of velocity
% nodes - see help for convertMeshToSecondOrder for detailed description of
% nVnodes) elements in the solution vector u correspond to x-velocity
% values in the subsequent nodes of the mesh, the next elements in u, this
% means u(nVnodes+1:2*nVnodes) contain values of the y-velocity for each
% mesh node, the following to z-velocity values and then the last elements,
% this means u(3*nVnodes+1:end) contain pressure values for each mesh node
% (note here that only first order nodes contain pressure solution - do not
% refer to mid-edge nodes when requesting pressure values, since this will
% cause an error with false indices in indices.indp vector). Using
% generatePressureData function might be helpful.
%
% Examples:
% 1. Extract y-velocity value for the 179th node in the mesh:
%       indices = generateIndices2D(p,t);
%       v179 = u(indices.indv(179));
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: CONVERTMESHTOSECONDORDER, EXPORTTOGMSH2D, GENERATEPRESSUREDATA,
%           IMPORTMESHGMSH.

% Pressure nodes
pNodesNumber = max(max(t(1:4, :)));

% Velocity nodes
vNodesNumber = size(p, 2);

% Save indices as logical vectors
indu = false(vNodesNumber*3 + pNodesNumber, 1);
indv = indu;
indw = indu;
indp = indu;

indu(1:vNodesNumber) = true;
indv(vNodesNumber + 1:2*vNodesNumber) = true;
indw(2*vNodesNumber + 1:3*vNodesNumber) = true;
indp((1:pNodesNumber) + 3*vNodesNumber) = true;

indices = struct('indp', indp, 'indu', indu, 'indv', indv, 'indw', indw);

end


%%



