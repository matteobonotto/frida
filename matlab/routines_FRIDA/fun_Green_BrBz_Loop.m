function [ Br,Bz ] = fun_Green_BrBz_Loop( source, point, curr) %#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Computation of Aphi, Br, Bz, Psi for axisymmetric loop current (no
%   thickness of the loop)
%
%   source (sources' geometry) - structure including:
%     - R: radial distance form the axis of the coil's centre [m]
%     - Z: vertical distance from z=0 plane of the coil's centre [m]
%
%   point (evaluation points) - structure including:
%     - RR: array of the radial coordinate of the evaluatin points [m]
%     - ZZ: array of the vertical coordinate of the evaluatin points [m]
%
%   res (results) - structure including:
%    - br (npt x ncoil)
%    - bz (npt x ncoil)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
npt=numel(point(:,1));
ncoil=numel(source(:,1));
br=zeros(npt,1);
bz=zeros(npt,1);

%% Computation (vectorization of quantities)

if ncoil <= npt
    
    RR=point(:,1);
    ZZ=point(:,2);
    for  jj=1:ncoil
        
        r0=source(jj,1);
        z0=source(jj,2);
        
        [res_br,res_bz] = fun_calcBrBz_1(r0,z0,RR,ZZ);
        
        br = br + res_br*curr(jj);
        bz = bz + res_bz*curr(jj);
        
    end
    
else
    
    r0=source(:,1);
    z0=source(:,2);
    for  jj=1:npt
        
        RR = point(jj,1);
        ZZ = point(jj,2);
        
        [res_br,res_bz] = fun_calcBrBz_2(r0,z0,RR,ZZ);
        
        tmp_r =res_br'*curr;
        tmp_z =res_bz'*curr;
        
        br(jj)=tmp_r;
        bz(jj)=tmp_z;
        
    end
end


%% Find points on axis (r=0)
RR=point(:,1);
ZZ=point(:,2);

ind_axis=find(point(:,1)==0);
br(ind_axis,:)=0;
if isempty(ind_axis)==0
    for  jj=1:ncoil
        r0=source(jj,1);
        z0=source(jj,2);
        bz(ind_axis,jj)=4.d-7*pi*r0/(2*(r0^2+(ZZ(ind_axis)-z0).^2).^1.5d0);
    end
end

%% Output
Br=br;
Bz=bz;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function  [res_br,res_bz] = fun_calcBrBz_1(r0,z0,RR,ZZ) %codegen

kk=sqrt(4*r0*RR./((r0+RR).^2+(ZZ-z0).^2));
kk_square = kk.^2;
kk_square(kk_square>1) = 1;

[J1,J2] = ellipke(kk_square);

res_br=+1.d-7.*kk.*(ZZ-z0)./(RR.*sqrt(r0*RR)).*...
    (-J1+(r0.^2+RR.^2+(ZZ-z0).^2)./((r0-RR).^2+(ZZ-z0).^2).*J2);
res_bz=+1.d-7.*kk./sqrt(r0*RR).*(J1+(r0^2-RR.^2-(ZZ-z0).^2)./...
    ((r0-RR).^2+(ZZ-z0).^2).*J2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function  [res_br,res_bz] = fun_calcBrBz_2(r0,z0,RR,ZZ) %codegen

kk=sqrt(4*r0*RR./((r0+RR).^2+(ZZ-z0).^2));
kk_square = kk.^2;
kk_square(kk_square>1) = 1;

[J1,J2] = ellipke(kk_square);

res_br=+1.d-7.*kk.*(ZZ-z0)./(RR*sqrt(r0*RR)).*...
    (-J1+(r0.^2+RR^2+(ZZ-z0).^2)./((r0-RR).^2+(ZZ-z0).^2).*J2);
res_bz=+1.d-7.*kk./sqrt(r0*RR).*(J1+(r0.^2-RR^2-(ZZ-z0).^2)./...
    ((r0-RR).^2+(ZZ-z0).^2).*J2);
end



















