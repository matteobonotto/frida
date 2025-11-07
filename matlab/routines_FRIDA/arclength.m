function [ARCLENGTH]=arclength(RR,ZZ)
%% [ARCLENGTH]=ARCLENGTH(X,Y)
%   LENGTH OF ARC DEFINED IN CARTESIAN COORDINATES
%       
%   M.B. 02/2018

%% CHEK IF THE CURVE IS OPEN OR CLOSEDrr=[RR; RR(1)];
rr=[RR; RR(1)];
zz=[ZZ; ZZ(1)];
distance=sqrt((rr(2:end,1)-rr(1:end-1,1)).^2+(zz(2:end,1)-zz(1:end-1,1)).^2);
dist_avg=sum(distance)./length(distance);
if RR(1)~=RR(end) && ZZ(1)~=ZZ(end)
    if max(distance) > 1.5*dist_avg
        % OPEN CURVE
        CLOSED = 0;
        [~,b]=max(distance);
        % %         figure; hold on; axis equal;
        % %         plot(RR(b),ZZ(b),'*');  plot(RR(b+1),ZZ(b+1),'*')
        index=[b+1:length(RR) 1:b]';
        RR=RR(index);
        ZZ=ZZ(index);
    else
        % CLOSED CURVE
        CLOSED = 1;
        RR=[RR; RR(1)];
        ZZ=[ZZ; ZZ(1)];
    end
end
npt=length(RR);

% % figure
% % plot(RR,ZZ,'.');
% % hold on
% % plot(RR(1), ZZ(1),'*')
% % plot(RR(end), ZZ(end),'*')

%% CUBIC SPLINE INTERPOLATION
Xspline=spline(1:npt,RR);
Yspline=spline(1:npt,ZZ);
Cx=Xspline.coefs;
Cy=Yspline.coefs;

tt=0:.05:1;
xx_t=zeros(npt-1,length(tt));
yy_t=zeros(npt-1,length(tt));
for ii=1:npt-1
    xx_t(ii,:)=Cx(ii,1)*tt.^3+Cx(ii,2)*tt.^2+Cx(ii,3)*tt+Cx(ii,4);
    yy_t(ii,:)=Cy(ii,1)*tt.^3+Cy(ii,2)*tt.^2+Cy(ii,3)*tt+Cy(ii,4);
end

%% ARC LENGTH VIA GAUSS QUADRATURE
[nodes,weight] = fun_GaussPoints_1D_MB(4);

nodes = nodes';
weight = weight';

a=0;
b=1;

arc=0;
tt_ii=.5*(b-a)*nodes+.5*(b+a);
for ii=1:npt-1
    xx_ii=3*Cx(ii,1)*tt_ii.^2+2*Cx(ii,2)*tt_ii+Cx(ii,3);
    yy_ii=3*Cy(ii,1)*tt_ii.^2+2*Cy(ii,2)*tt_ii+Cy(ii,3);
    
    ff_ii=(xx_ii.^2+yy_ii.^2).^.5;
    arc_ii=.5*(b-a)*ff_ii*weight';
    arc=arc+arc_ii;
end

%% PROVA
distance=sqrt((RR(2:end,1)-RR(1:end-1,1)).^2+(ZZ(2:end,1)-ZZ(1:end-1,1)).^2);
dist=sum(distance);
% % 100*(arc-dist)/arc
if 100*(arc-dist)/arc > 1
    error(['FAILED INTEGRATION, RELATIVE ERROR = ' num2str(100*(arc-dist)/arc)])
else
    ARCLENGTH=arc;
end

end

