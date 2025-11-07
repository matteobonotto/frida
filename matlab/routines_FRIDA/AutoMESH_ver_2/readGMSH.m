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


