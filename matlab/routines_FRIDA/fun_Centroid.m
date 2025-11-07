function [solk1]=fun_Centroid(meshData,solk1)

ind_n_INpla=meshData.ind_n_INpla;
% Zc
Zc=sum(meshData.n(ind_n_INpla,2).*solk1.Iphi(ind_n_INpla))/...
        sum(solk1.Iphi(ind_n_INpla));

% Rc
Rc=sqrt(sum(meshData.n(ind_n_INpla,1).^2.*solk1.Iphi(ind_n_INpla))/...
    sum(solk1.Iphi(ind_n_INpla)));

distanza=sqrt((Rc-meshData.n(ind_n_INpla,1)).^2+...
    (Zc-meshData.n(ind_n_INpla,2)).^2);
[~,ind_min]=min(distanza);

solk1.index_centroid=ind_n_INpla(ind_min);
% % solk1.Centroid_ZZ=meshData.n(ind_n_INpla(ind_min),2);
% % solk1.Centroid_RR=meshData.n(ind_n_INpla(ind_min),1);
solk1.Centroid_ZZ = Zc ;
solk1.Centroid_RR = Rc;

%%%
% % ind_t_Centroid = 
% % psi_C = interp1()

% % plot(Rc,Zc,'*r');
% % plot(solk1.Centroid_RR,solk1.Centroid_ZZ,'ob');


distanza=sqrt((Rc-meshData.nodes_Gauss_pla(:,1)).^2+...
    (Zc-meshData.nodes_Gauss_pla(:,2)).^2);
[~,ind_min]=min(distanza);

ind_t_Centroid = ceil(ind_min/meshData.n_Gauss);
ind_n_centroid = meshData.t(ind_t_Centroid,:);

N_order = meshData.shape_functions_N_order;
cc_tri = reshape(meshData.shape_functions(ind_t_Centroid,:),3*N_order,3*N_order).';

Psi_Centroid = solk1.Psi(ind_n_centroid).'*(cc_tri*([Rc^2 Zc^2 Rc*Zc Rc Zc 1].'));
solk1.Psi_Centroid = Psi_Centroid;
% % figure
% % plot3(meshData.n(ind_n_centroid,1),meshData.n(ind_n_centroid,2),solk1.Psi(ind_n_centroid),'o')
% % hold on;
% % plot3(Rc,Zc,Psi_Centroid,'*')

