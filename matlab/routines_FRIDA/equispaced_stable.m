function [RR_EQUI,ZZ_EQUI]=equispaced_stable(RR,ZZ,factor)

if norm([RR(1) ZZ(1)] - [RR(end) ZZ(end)]) < eps
    CLOSEDCURVE = true;
else
    CLOSEDCURVE = false;
    RR = [RR; RR(1)];
    ZZ = [ZZ; ZZ(1)];
end

%%
NPT=length(RR);
RR=interp1(RR,linspace(1,NPT,10000))';
ZZ=interp1(ZZ,linspace(1,NPT,10000))';


NPT_EQUI=length(RR);

DIST=sqrt((RR(2:end)-RR(1:end-1)).^2+(ZZ(2:end)-ZZ(1:end-1)).^2);
LENGTH = sum(DIST);


s_coo=zeros(NPT_EQUI,1);
for ii=2:NPT_EQUI
    s_coo(ii)=s_coo(ii-1)+DIST(ii-1);
end

DIST_NEW=LENGTH/factor;

RR_EQUI=interp1(s_coo,RR(1:end),(s_coo(1):DIST_NEW:s_coo(end)))';
ZZ_EQUI=interp1(s_coo,ZZ(1:end),(s_coo(1):DIST_NEW:s_coo(end)))';

%%
if CLOSEDCURVE == false
    RR_EQUI = RR_EQUI(1:end-1);
    ZZ_EQUI = ZZ_EQUI(1:end-1);
end












