function [lati,lati_bordo]=fun_lati(meshData,LATIBORDO)
% find indices of edges, no duplicates
% find boundary edges

%% Edges
tri_srt=meshData.t_choosen;
ntri=size(meshData.t_choosen,1);
lati=zeros(1,2);
lati_tot=zeros(3*size(meshData.t_choosen,1),2);
for i=1:ntri
    triangle=tri_srt(i,:);
    lati_loc=[triangle(1) triangle(2);
        triangle(2) triangle(3);
        triangle(3) triangle(1)];
    if i==1
        lati=lati(2:end,:);
    end
    lati_tot(3*i-2:3*i,:)=lati_loc;
end
lati_tot=sort(lati_tot,2);
[lati,~,~]=unique(lati_tot,'rows');

%% Boundary edges
lati_bordo=zeros(size(lati));
if LATIBORDO==1
    aa=sort(meshData.t_choosen(:,[1 2]),2);
    bb=sort(meshData.t_choosen(:,[2 3]),2);
    cc=sort(meshData.t_choosen(:,[3 1]),2);
    tic
    parfor i=1:size(lati,1)
        lati_i=lati(i,:);
        ind1=find(ismember(aa,lati_i,'rows')==1);
        ind2=find(ismember(bb,lati_i,'rows')==1);
        ind3=find(ismember(cc,lati_i,'rows')==1);
        ind=[ind1; ind2; ind3];
        if numel(ind)==1
            lati_bordo(i,:)=lati(i,:);
            % %         lati_bordo(ii,:)=lati(i,:);
            % %         ii=ii+1;
        end
    end
    toc
    lati_bordo=lati_bordo(sum(lati_bordo,2)~=0,:);
end
end

