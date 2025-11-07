function [shape_f,shape_f_norm,error_inv,cond_AA] = fun_shape_functions_stable(tri,nodes,N_order) %#codegen


nt = size(tri,1);

shape_f = zeros(nt,3*N_order*3*N_order);
error_1 = zeros(nt,1);
error_2 = zeros(nt,1);

error_inv = zeros(nt,1);
cond_AA = zeros(nt,1);
for ii = 1:nt
    tri_ii = tri(ii,:);
    
    if N_order == 1
        
        P1 = nodes(tri_ii(:,1),:);
        P2 = nodes(tri_ii(:,2),:);
        P3 = nodes(tri_ii(:,3),:);
        
        AA = [P1 1; P2 1; P3 1];
        
        coeff = inv(AA);
        shape_f(ii,:) = reshape(coeff,1,(3*N_order)^2);
        
        error_inv(ii) = norm(coeff*AA-eye(3));
        
        cond_AA(ii) = cond(AA);
                
    elseif  N_order == 2
        
        P1 = nodes(tri_ii(:,1),:);
        P2 = nodes(tri_ii(:,2),:);
        P3 = nodes(tri_ii(:,3),:);
        P4 = nodes(tri_ii(:,4),:);
        P5 = nodes(tri_ii(:,5),:);
        P6 = nodes(tri_ii(:,6),:);
        
            AA = [P1.^2 P1(1)*P1(2) P1 1; ...
                P2.^2 P2(1)*P2(2) P2 1; ...
                P3.^2 P3(1)*P3(2) P3 1; ...
                P4.^2 P4(1)*P4(2) P4 1; ...
                P5.^2 P5(1)*P5(2) P5 1; ...
                P6.^2 P6(1)*P6(2) P6 1];
        
        
        coeff_1 = (AA)\eye(size(AA));
        coeff_2 = inv(AA);
        
        error_1(ii) = norm(coeff_1*AA-eye(6));
        error_2(ii) = norm(coeff_2*AA-eye(6));
        
        error_inv(ii) = min([error_1(ii) error_2(ii)]);
        
        cond_AA(ii) = cond(AA);
        
        if norm(coeff_1*AA-eye(6)) <= norm(coeff_2*AA-eye(6))
            shape_f(ii,:) = reshape(coeff_1,1,(3*N_order)^2);
        else
            shape_f(ii,:) = reshape(coeff_2,1,(3*N_order)^2);
        end        
        
    end
    
end


%%

shape_f_norm = [];

if N_order == 1
    
    P1 = [0 0];
    P2 = [1 0];
    P3 = [0 1];
    
    AA_norm = [P1 1; P2 1; P3 1];
    
    coeff_norm = inv(AA_norm);
    
    shape_f_norm = reshape(coeff_norm,1,(3*N_order)^2);
    
elseif  N_order == 2
    
    P1 = [0 0];
    P2 = [1 0];
    P3 = [0 1];
    P4 = [0.5 0];
    P5 = [.5 .5];
    P6 = [0 .5];
    
    AA_norm = [P1.^2 P1(1)*P1(2) P1 1; ...
        P2.^2 P2(1)*P2(2) P2 1; ...
        P3.^2 P3(1)*P3(2) P3 1; ...
        P4.^2 P4(1)*P4(2) P4 1; ...
        P5.^2 P5(1)*P5(2) P5 1; ...
        P6.^2 P6(1)*P6(2) P6 1];
    
    coeff_norm = inv(AA_norm);
    
    shape_f_norm = reshape(coeff_norm,1,(3*N_order)^2);
    
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




















