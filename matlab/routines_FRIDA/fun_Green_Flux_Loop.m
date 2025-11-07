function [ psi ] = fun_Green_Flux_Loop( source, point, curr) %#codegen
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
%    - a (npt x ncoil)
%    - br (npt x ncoil)
%    - bz (npt x ncoil)
%    - psi (npt x ncoil)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

npt   = numel(point(:,1));
ncoil = numel(source(:,1));
psi   = zeros(npt,1);

% % psi=zeros(npt,ncoil);

%% Computation (vectorization of quantities)

if ncoil <= npt
    
    RR=point(:,1);
    ZZ=point(:,2);
    
    for  jj=1:ncoil
        
        r0=source(jj,1);
        z0=source(jj,2);
        
        kk=sqrt(4*r0*RR./((r0+RR).^2+(ZZ-z0).^2));
        kk_square = kk.^2;
        kk_square(kk_square>1) = 1;
        
        [J1,J2] = ellipke(kk_square);
        
        res_psi=4.d-7./kk.*sqrt(r0./RR).*((1-kk.^2/2).*J1-J2).*(2*pi*RR);
        
        psi = psi + res_psi*curr(jj);
        
    end
    
else
    
    r0=source(:,1);
    z0=source(:,2);
    
    for  jj=1:npt
        
        RR = point(jj,1);
        ZZ = point(jj,2);
    
        kk=sqrt(4*r0*RR./((r0+RR).^2+(ZZ-z0).^2));
        kk_square = kk.^2;
        kk_square(kk_square>1) = 1;
        
        [J1,J2] = ellipke(kk_square);
        
        res_psi=4.d-7./kk.*sqrt(r0./RR).*((1-kk.^2/2).*J1-J2).*(2*pi*RR);
        
        tmp_psi =res_psi'*curr;
        
        psi(jj)=tmp_psi;
        
    end
end

%% Find points on axis (r=0)
psi(point(:,1)==0,:)=0;

%% Output
% % res.psi=psi;

end

