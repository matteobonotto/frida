function [lati,lati_bordo]=fun_lati_fast(t_choosen,LATIBORDO) %codegen
% find indices of edges, no duplicates
% find boundary edges

%% Edges
tri_srt=t_choosen;
ntri=size(t_choosen,1);
npoints = size(tri_srt,2);
lati_tot=zeros(npoints*size(t_choosen,1),2);
for i=1:ntri
    triangle=tri_srt(i,:);
    ind_lati_loc = [1 2; 2 3; 3 1];
    
    lati_loc=triangle(ind_lati_loc);
    
    lati_tot(3*(i-1)+1:3*i,:)=lati_loc;
end
lati_tot=sort(lati_tot,2);
[lati,~,~]=unique(lati_tot,'rows');

%% Boundary edges
lati_bordo=zeros(size(lati));
if LATIBORDO==1
    
    % %     if npoints == 3
    
    aa=sort(t_choosen(:,[1 2]),2);
    bb=sort(t_choosen(:,[2 3]),2);
    cc=sort(t_choosen(:,[3 1]),2);
    % %         tic
    for i=1:size(lati,1) % NOTE: parfor only in the mex version
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
    % %         toc
    lati_bordo=lati_bordo(sum(lati_bordo,2)~=0,:);
    
    
    
    
end




end

