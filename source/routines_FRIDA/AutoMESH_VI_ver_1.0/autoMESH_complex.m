function [p,e,t]=autoMESH_complex(filename,GEOMETRY,order,dirname)

pp = pwd;
% % cd(dirname)


%% WRITEB INPUT GMSH .geo FILE
writeGMSH(filename,GEOMETRY);


%% MESH THE GEOMETRY USING GMSH;

command_dos = sprintf('gmsh -2 -order %i -format m %s.geo', order, filename);

if ismac
    setenv('PATH', getenv('PATH')+":/opt/homebrew/bin/")
    message = ['autoMESH: if using MAC add gmsh dir ', ...
        'like this -> ', ...
        'setenv("PATH", getenv("PATH")+":/opt/homebrew/bin/")'];
    warning(message)
end

dos(command_dos);

%% READ MESH

run(filename)

p = msh.POS(:,[1 2]);
if order == 1
    t = msh.TRIANGLES(:,1:end-1);
elseif order == 2
    t = msh.TRIANGLES6(:,1:end-1);
elseif order == 3
    t = msh.TRIANGLES10(:,1:end-1);
end

% EDGES
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


cd(pp)


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
end


n_surf  = 1;
for ii=1:size(surface,2)

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




