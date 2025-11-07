function [Er,Ar,Br,Cr,Dr] = funBlockKrylov_out(E,A,B,C,D,q2,s2)
%KRYLOV_MOR
% Krylov per problemi MIMO, twosided-Arnoldi

BB = C';
nn = size(A,2);
mm = size(BB,2);
[LAA,UAA] = lu(-A+s2*E);
invAA = inv(UAA)*inv(LAA);
AW1 = invAA';
BW1 = invAA'*BB;

% onesided-Arnoldi
[W] = funArnoldi(AW1,BW1,q2);
V = W;

% Reduction
Er = W'*E*V;
Ar = W'*A*V;
Br = W'*B;
Cr = C*V;
Dr = D;

end

