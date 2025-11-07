function [MatInd_nodes_tri] = fun_MatInd_nodes_tri(N_order,tri,nn) %#codegen

MatInd_nodes_tri = [];

if N_order == 1
    
% %     num_el = zeros(nn,1);
% %     for ii=1:nn
% %         qq = unique([find(tri(:,1) == ii); ...
% %             find(tri(:,2) == ii); ...
% %             find(tri(:,3) == ii)]);
% %         
% %         num_el(ii) = numel(qq);
% %         
% %     end
% %     
% %     num_el_max = max(num_el);
    
    MatInd_nodes_tri = zeros(nn,5);
    for ii=1:nn
        qq = unique([find(tri(:,1) == ii); ...
            find(tri(:,2) == ii); ...
            find(tri(:,3) == ii)]);
        
    MatInd_nodes_tri(ii,1:numel(qq)) = qq;
        
    end
    
elseif N_order == 2
    
% %     num_el = zeros(nn,1);
% %     for ii=1:nn
% %         qq = unique([find(tri(:,1) == ii); ...
% %             find(tri(:,2) == ii); ...
% %             find(tri(:,3) == ii); ...
% %             find(tri(:,4) == ii); ...
% %             find(tri(:,5) == ii); ...
% %             find(tri(:,6) == ii)]);
% %         
% %         num_el(ii) = numel(qq);
% %         
% %     end
% %     
% %     num_el_max = max(num_el);
    
    MatInd_nodes_tri = zeros(nn,8);
    for ii=1:nn
        qq = unique([find(tri(:,1) == ii); ...
            find(tri(:,2) == ii); ...
            find(tri(:,3) == ii); ...
            find(tri(:,4) == ii); ...
            find(tri(:,5) == ii); ...
            find(tri(:,6) == ii)]);
        
    MatInd_nodes_tri(ii,1:numel(qq)) = qq;
        
    end
    
end

