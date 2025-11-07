function  [Area_coils]=fun_AreaCoils_FEM(meshData,Conductors)

% % figure; hold on
Area_coils = zeros(Conductors.Nconductors,1);
for ii = 1:Conductors.Nconductors
    
    tri_ii = meshData.t(meshData.type == ii,1:3);
    
% %     triplot(meshData.t(meshData.type == ii,1:3),meshData.n(:,1),meshData.n(:,2))
    
    P1 = meshData.n(tri_ii(:,1),:);
    P2 = meshData.n(tri_ii(:,2),:);
    P3 = meshData.n(tri_ii(:,3),:);
    
    edge_1 = P2 - P1;
    edge_2 = P3 - P2;

    Area_ii = .5*abs(edge_1(:,1).*edge_2(:,2) - edge_1(:,2).*edge_2(:,1));
    
    Area_coils(ii) = sum(Area_ii);
% %     pause
end




