function [shape_f_norm,nodes_norm,N_nodes_sf] = fun_shape_functions_norm(N_order) %#codegen



%%

shape_f_norm = [];

if N_order == 1
        
    nodes_norm = [0 0;
        1 0;
        0 1];
    
    % %     AA_norm = [nodes_norm ...
    % %         ones(3*N_order,1)];
    % %
    % %     coeff_norm = inv(AA_norm);
    
    coeff_norm = [ ...
    -1     1     0
    -1     0     1
     1     0     0];
    
    shape_f_norm = reshape(coeff_norm,1,(3*N_order)^2);
    
elseif  N_order == 2
    
    nodes_norm = [0 0;
        1 0;
        0 1;
        0.5 0;
        .5 .5;
        0 .5];
    
    % %     AA_norm = [nodes_norm.^2 ...
    % %         nodes_norm(:,1).*nodes_norm(:,2) ...
    % %         nodes_norm ...
    % %         ones(3*N_order,1)];
    
    % %     coeff_norm = inv(AA_norm);
    
    coeff_norm = [ ...
     2     2     0    -4     0     0
     2     0     2     0     0    -4
     4     0     0    -4     4    -4
    -3    -1     0     4     0     0
    -3     0    -1     0     0     4
     1     0     0     0     0     0];
    
 % %     norm(coeff_norm*AA_norm   eye(size(AA_norm)))
    
    shape_f_norm = reshape(coeff_norm,1,(3*N_order)^2);
    
elseif  N_order == 3
    
    nodes_norm = [0 0;
        1 0;
        0 1;
        1/3 0;
        2/3 0;
        2/3 1/3;
        1/3 2/3;
        0 2/3;
        0 1/3;
        1/3 1/3];
        
    % %     AA_norm = [nodes_norm.^3 ...
    % %         nodes_norm(:,1).^2.*nodes_norm(:,2) ...
    % %         nodes_norm(:,1).*nodes_norm(:,2).^2 ...
    % %         nodes_norm.^2 ...
    % %         nodes_norm(:,1).*nodes_norm(:,2) ...
    % %         nodes_norm ...
    % %         ones(size(nodes_norm,1),1)];
    % %
    % %
    % %     coeff_norm = inv(AA_norm);
    
      coeff_norm = [...
       -4.5000    4.5000         0   13.5000  -13.5000         0         0         0         0         0
       -4.5000         0    4.5000         0         0         0         0  -13.5000   13.5000         0
      -13.5000         0         0   27.0000  -13.5000   13.5000         0         0   13.5000  -27.0000
      -13.5000         0         0   13.5000         0         0   13.5000  -13.5000   27.0000  -27.0000
        9.0000   -4.5000         0  -22.5000   18.0000         0         0         0         0         0
        9.0000         0   -4.5000         0         0         0         0   18.0000  -22.5000         0
       18.0000         0         0  -22.5000    4.5000   -4.5000   -4.5000    4.5000  -22.5000   27.0000
       -5.5000    1.0000         0    9.0000   -4.5000         0         0         0         0         0
       -5.5000         0    1.0000         0         0         0         0   -4.5000    9.0000         0
        1.0000         0         0         0         0         0         0         0         0         0];
    
% %     norm(coeff_norm*AA_norm   eye(size(AA_norm)))


    shape_f_norm = reshape(coeff_norm,1,(3*N_order)^2);
    
else
    
    shape_f_norm = NaN;
    nodes_norm= NaN;
    N_nodes_sf = NaN;
    
end


N_nodes_sf = size(nodes_norm,1);





end




% % P1 = nodes(tri(:,1),:);
% % P2 = nodes(tri(:,2),:;
% % P3 = nodes(tri(:,3),:);
% % center = [sum([P1(:,1) P2(:,1) P3(:,1)],2) sum([P1(:,2) P2(:,2) P3(:,2)],2)]/3;
% %
% %
% % tic
% % F1 = scatteredInterpolant(center(:,1),center(:,2),error_1);
% % F2 = scatteredInterpolant(center(:,1),center(:,2),error_2);
% % F3 = scatteredInterpolant(center(:,1),center(:,2),error_3);
% % J_interp1 = F1(nodes(:,1),nodes(:,2));
% % J_interp2 = F2(nodes(:,1),nodes(:,2));
% % J_interp3 = F3(nodes(:,1),nodes(:,2));
% % toc
% %
% %
% %
% %
% % Jpla1.faces=tri(:,1:3);
% % Jpla1.vertices=nodes;
% % Jpla1.facevertexcdata=J_interp1;
% % Jpla2.faces=tri(:,1:3);
% % Jpla2.vertices=nodes;
% % Jpla2.facevertexcdata=J_interp2;
% % Jpla3.faces=tri(:,1:3);
% % Jpla3.vertices=nodes;
% % Jpla3.facevertexcdata=J_interp3;
% % Jpla4.faces=tri(:,1:3);
% % Jpla4.vertices=nodes;
% % qq = zeros(size(J_interp1));
% % qq(J_interp1<J_interp2) = 1;
% % Jpla4.facevertexcdata=qq;
% %
% %
% % figure
% % subplot(1,3,1)
% % hh=patch(Jpla1,'facecolor','interp','edgecolor','none');
% % caxis([min([error_1; error_2; error_3]) max([error_1; error_2; error_3])])
% % axis square; box on; colorbar vert; colormap jet
% % title('lu(A)')
% % subplot(1,3,2)
% % hh=patch(Jpla2,'facecolor','interp','edgecolor','none');
% % caxis([min([error_1; error_2; error_3]) max([error_1; error_2; error_3])])
% % title('inv(A)')
% % axis square; box on; colorbar vert; colormap jet
% % subplot(1,3,3)
% % hh=patch(Jpla3,'facecolor','interp','edgecolor','none');
% % caxis([min([error_1; error_2; error_3]) max([error_1; error_2; error_3])])
% % title('best')
% % axis square; box on; colorbar vert; colormap jet
% %
% %
% % figure
% % hh=patch(Jpla4,'facecolor','interp','edgecolor','none');
% % title('best')
% % axis square; box on; colorbar vert; colormap jet




















