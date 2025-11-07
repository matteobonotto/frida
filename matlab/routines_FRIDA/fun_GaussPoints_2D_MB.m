function [xx,yy,ww] = fun_GaussPoints_2D_MB(nn)

%% Define Gauss quadrature rule for a rectangle (-1,1)x(-1,1)
[xx,ww] = fun_GaussPoints_1D_MB(nn);


% compute Gauss quadrature rule for 2D
[XX,YY] = meshgrid(xx,xx); 
WW = ww*ww';

xx = XX(:);  
yy = YY(:);  
ww = WW(:);



























