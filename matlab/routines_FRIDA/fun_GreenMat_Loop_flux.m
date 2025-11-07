function [ res ] = fun_GreenMat_Loop_flux( source, point ) %#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Computation of Psi for axisymmetric loop current (no
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
%    - psi (npt x ncoil)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
npt=numel(point.RR);
ncoil=numel(source.R); 
% % a=zeros(npt,ncoil);
% % br=zeros(npt,ncoil);
% % bz=zeros(npt,ncoil);
psi=zeros(npt,ncoil);
RR=point.RR;
ZZ=point.ZZ;

%% Computation (vectorization of quantities)
for  jj=1:ncoil
    r0=source.R(jj);
    z0=source.Z(jj);
    kk=sqrt(4*r0*RR./((r0+RR).^2+(ZZ-z0).^2));
    kk_square = kk.^2;
    kk_square(kk_square>1) = 1;
    [J1,J2] = ellipke(kk_square);
    res_psi=4.d-7./kk.*sqrt(r0./RR).*((1-kk.^2/2).*J1-J2).*(2*pi*RR);
    psi(:,jj)=res_psi;
end

%% Find points on axis (r=0)
% % ind_axis=find(point.RR==0);
psi(point.RR==0,:)=0;

%% Output

res.psi=psi;

end

