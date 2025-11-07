function [Er,Ar,Br,Cr,Dr,V] = funBlockKrylov_in(E,A,B,C,D,q1,s1)
%KRYLOV_MOR
% Krylov per problemi MIMO, twosided-Arnoldi

nn = size(A,2);
mm = size(B,2);

I = eye(size(A,1));
[LA,UA] = lu((-A+s1*E));
invA = inv(UA)*inv(LA);
invAB = invA*B;

AV1 = invA;
BV1 = invAB;

% twosided-Arnoldi
[V] = funArnoldi(AV1,BV1,q1);
% % [V] = funArnoldi_refined(AV1,BV1,q1);
W = V;

% % [W] = funArnoldi(A',C',q,ss);
% % V = W;

% Reduction
Er = W'*E*V;
Ar = W'*A*V;
Br = W'*B;
Cr = C*V;

if ~isempty(D)
    Dr = D;
end


end

