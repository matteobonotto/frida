function [p,e,t]=autoMESH(filename,GEOMETRY,SPACING)

%% WRITEB INPUT GMSH .geo FILE 
writeGMSH(filename,GEOMETRY,SPACING);

%% MESH THE GEOMETRY USING GMSH
dos(['gmsh -2 ' filename '.geo']);

%% READ MESH
[p,e,t]=readGMSH(filename);

end


%%
function [ fileID ] = writeGMSH(filename,GEOMETRY,SPACING)

if abs(GEOMETRY(1,1)-GEOMETRY(end,1)) < 1E-10 && ...
     abs(GEOMETRY(1,2)-GEOMETRY(end,2)) < 1E-10
    GEOMETRY=GEOMETRY(1:end-1,:);
end

fileID = fopen([filename '.geo'],'w');

shape{1}=GEOMETRY;
element=shape{1};
spacing=zeros(size(element,1),1);
spacing(:)=SPACING;

%% POINTS
n_point=1;
for ii = 1:size(element,1)
    fprintf(fileID,'%s  %d %s %8f %s %8f %s %8f %s \n',...
        char('Point('),...
        n_point,...
        char(')={'),...
        element(ii,1),...
        ',',...
        element(ii,2),...
        ',0,',...
        spacing(ii),...
        '};');
    n_point=n_point+1;
end

fprintf(fileID,'%6s \n',' ');

%% LINES
n_line=1;
for ii = 1:size(element,1)
    if ii<size(element,1)
        fprintf(fileID,'%s  %d %s %d %s %d %s \n',...
            char('Line('),...
            n_line,...
            char(')={'),...
            n_line,...
            ',',...
            n_line+1,...
            '};');
    elseif ii==size(element,1)
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

%% LINE LOOP
n_loop=1;
fprintf(fileID,'%s \n','Line Loop(1)={');
for ii = 1:size(element,1)
    if ii<size(element,1)
        fprintf(fileID,'%d %s\n',...
            n_line-size(element,1)+ii-1,...
            ',');
    elseif ii==size(element,1)
        fprintf(fileID,'%d %s\n',...
            n_line-size(element,1)+ii-1, ...
            '};');
    end
end

fprintf(fileID,'%1s \n',' ');

%% PLANE SURFACE
fprintf(fileID,'%s \n','Plane Surface(1)=(1);');


fclose(fileID);

end

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

