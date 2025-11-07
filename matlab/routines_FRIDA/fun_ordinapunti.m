function [gamma_sort,ind_TH] = fun_ordinapunti(gamma,pole)
%FUN_ORDINAPUNTI sort a given set of point in a clockwise order

%%
if nargin == 1
    baricentro=[sum(gamma(:,1))  sum(gamma(:,2))]/size(gamma,1);
else
    baricentro = pole;
end
gamma_centro(:,1)=gamma(:,1)-baricentro(1);
gamma_centro(:,2)=gamma(:,2)-baricentro(2);

[TH,~]=cart2pol(gamma_centro(:,1),gamma_centro(:,2));
TH(TH < 0) = TH(TH < 0)+2*pi;
[~,ind_TH]=sort(TH);

gamma_sort=gamma(ind_TH,:);

%%
% % figure
% % hold on
% % plot(gamma(:,1),gamma(:,2),'.')
% % plot(baricentro(1),baricentro(2),'*')
% % for i=1:size(gamma,1)
% %     plot(gamma_sort(i,1),gamma_sort(i,2),'*')
% %     pause
% % end

end

