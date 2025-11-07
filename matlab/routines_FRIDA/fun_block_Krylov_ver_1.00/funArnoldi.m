function [VV] = funArnoldi(K1,K2,q)
% Block-oneside-ARNOLDI for MIMO systems
% simple selection procedure

% rimozione vettori linearmente dipendenti
if size(K2,2) > 1    
    B_m1 = licols(K2,1e-12);
else
    B_m1 = K2;
end
m1 = size(B_m1,2);
VV = [];

% first starting vector
VV(:,1) = B_m1(:,1)/sqrt(B_m1(:,1)'*B_m1(:,1));

for i = 2:q
    clear r v_hat h v norm_v_hat;
    if i <= m1
       v_hat = B_m1(:,i);
    else
       v_hat = K1*VV(:,i-m1);
    end    
    % orthogonalization
    for j = 1:i-1
        h = VV(:,j)'*v_hat;
        v_hat = v_hat - h*VV(:,j);
    end 
    % check if vector is zero
% %     norm_v_hat = norm(v_hat);
    % normalization
    v = v_hat/sqrt(v_hat'*v_hat);
    VV(:,i) = v;
end


end



































% % % 
% % % 
% % % 
% % % 
% % % Kq1 = invA*B;
% % % Kq2 = invA'*C';
% % % Kq1 = licols(Kq1,1e-12);
% % % Kq2 = licols(Kq2,1e-12);
% % % mm1 = size(Kq1,2);
% % % mm2 = size(Kq2,2);
% % % V = [];
% % % W = [];
% % % 
% % % % INPUT Krylov subspace
% % % V(:,1) = Kq1(:,1);
% % % V(:,1) = Kq1(:,1)/sqrt(Kq1(:,1)'*Kq1(:,1));
% % % rankV = rank(V);
% % % 
% % % i = 2;
% % % while rankV == size(V,2) %& i<=q 
% % % 
% % %     if i <= mm1
% % %         R = Kq1(:,i);   
% % %         v_hat = R;
% % %     else
% % %         clear R v_i v_hat val index h d
% % %         R = invA*V(:,i-mm1);
% % %         v_hat = R;
% % %     end
% % %     
% % %     for j = 1:i-1
% % %         clear h
% % %         h = v_hat'*V(:,j);
% % %         v_hat = v_hat - h*V(:,j);
% % %     end
% % %     v_i = v_hat;
% % %     V(:,i) = v_i/norm(v_i);
% % %     rankV = rank(V);
% % %     i = i+1;
% % % end
% % % 
% % % qv_best = size(V,2) -1;
% % % 
% % % 
% OUTPUT Krylov subspace
% % % W(:,1) = Kq2(:,1);
% % % W(:,1) = Kq2(:,1)/sqrt(Kq2(:,1)'*Kq2(:,1));
% % % rankW = rank(W);
% % % 
% % % i = 2;
% % % while rankW == size(W,2) %& i<=q 
% % %     
% % %     if i <= mm2
% % %         R = Kq2(:,i);   
% % %         w_hat = R;       
% % %     else
% % %         clear R v_i v_hat val index h d
% % %         R = invA'*W(:,i-mm2);
% % %         w_hat = R;
% % %     end
% % %     
% % %     for j = 1:i-1
% % %         clear h
% % %         h = w_hat'*W(:,j);
% % %         w_hat = w_hat - h*W(:,j);
% % %     end
% % %     w_i = w_hat;
% % %     W(:,i) = w_i/norm(w_i);
% % %     rankW = rank(W);
% % %     i = i+1;
% % % end
% % % 
% % % 
% % % qw_best = size(W,2) -1;
% % % 
% % % disp 'Best reducing factor q: '
% % % qbest = min([qv_best qw_best])
% % % V = V(:,1:qbest);
% % % W = W(:,1:qbest);
% % % 
