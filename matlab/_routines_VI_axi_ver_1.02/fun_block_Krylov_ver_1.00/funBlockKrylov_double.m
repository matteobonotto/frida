function [Er,Ar,Br,Cr,Dr] = funBlockKrylov_double(E,A,B,C,D,qr_in,qr_out,s_in,s_out)
%KRYLOV_MOR
% Krylov per problemi MIMO, twosided-Arnoldi

nn = size(A,2);
mm = size(B,2);

%% INPUT
As1E = (-A+s_in*E);
AV1 = As1E\E;
BV1 = As1E\B;

[V] = funArnoldi(AV1,BV1,qr_in);

%% OUTPUT
AsTE = (-A+s_out*E)';
AV1 = AsTE\E';
BV1 = AsTE\C';

% onesided-Arnoldi
[W] = funArnoldi(AV1,BV1,qr_out);

%% Reduction
Er = W'*E*V;
Ar = W'*A*V;
Br = W'*B;
Cr = C*V;
Dr = D;

end

