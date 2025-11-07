function [Area_t,Area_n]=fun_AreaMesh_FEM(meshData)


N_order = meshData.shape_functions_N_order;

%% Area Trianglles

Area_t=zeros(size(meshData.t,1),1);

for i=1:size(meshData.t,1)
    tri=meshData.t(i,:);
    
    if N_order == 1
        P1=meshData.n(tri(1),:);
        P2=meshData.n(tri(2),:);
        P3=meshData.n(tri(3),:);
        
        triangle = [P1; P2; P3];
        
    elseif N_order == 2
        P1=meshData.n(tri(1),:);
        P2=meshData.n(tri(2),:);
        P3=meshData.n(tri(3),:);
        P4=meshData.n(tri(4),:);
        P5=meshData.n(tri(5),:);
        P6=meshData.n(tri(6),:);
        
        triangle = [P1; P2; P3; P4; P5; P6];
        [triangle] = fun_ordinapunti(triangle);
        
    end
    
    
    % %     figure
    % %     plot(triangle(:,1),triangle(:,2),'-r')
    
    % %     100*(Area_p(i) - polyarea(triangle(:,1),triangle(:,2)))/polyarea(triangle(:,1),triangle(:,2))
    Area_t(i)=polyarea(triangle(:,1),triangle(:,2));
    
end

%% Area nodes

if nargout > 1
    N_order = meshData.shape_functions_N_order;
    
    ww_Gauss = meshData.ww_Gauss_pla;
    P_Gauss = meshData.nodes_Gauss_pla;
    
    ind_t_selected = 1:meshData.nt;
    triangles = meshData.t(ind_t_selected,:);
    
    
    shape_functions = meshData.shape_functions;
    
    % % tic
    Area_n = zeros(meshData.nn,1);
    
    if N_order == 1
        
        for ii = 1:size(triangles,1)
            
            ind_t_ii = ind_t_selected(ii);
            tri_ii = triangles(ii,:);
            
            ii_tri_Gauss = (ind_t_ii-1)*3*N_order+1:ind_t_ii*3*N_order;
            
            ww = ww_Gauss(ii_tri_Gauss,:);
            nodes = P_Gauss(ii_tri_Gauss,:);
            
            cc = reshape(shape_functions(ind_t_ii,:),3*N_order,3*N_order)';
            
            tmp = [nodes ones(3*N_order,1)]';
            
            NN_1 = cc(1,:)*tmp;
            NN_2 = cc(2,:)*tmp;
            NN_3 = cc(3,:)*tmp;
             
            A_tri_1 = ww'*NN_1';
            A_tri_2 = ww'*NN_2';
            A_tri_3 = ww'*NN_3';
            
            A_tri = [A_tri_1; A_tri_2; A_tri_3];
            
            Area_n(tri_ii) = Area_n(tri_ii) + A_tri;
            
        end
        
    elseif N_order == 2
        
        for ii = 1:size(triangles,1)
            
            ind_t_ii = ind_t_selected(ii);
            tri_ii = triangles(ii,:);
            
            ii_tri_Gauss = (ind_t_ii-1)*3*N_order+1:ind_t_ii*3*N_order;
            
            ww = ww_Gauss(ii_tri_Gauss,:);
            nodes = P_Gauss(ii_tri_Gauss,:);
            
            cc = reshape(shape_functions(ind_t_ii,:),3*N_order,3*N_order)';
            
            tmp = [nodes.^2 nodes(:,1).*nodes(:,2) nodes ones(3*N_order,1)]';
            
            NN_1 = cc(1,:)*tmp;
            NN_2 = cc(2,:)*tmp;
            NN_3 = cc(3,:)*tmp;
            NN_4 = cc(4,:)*tmp;
            NN_5 = cc(5,:)*tmp;
            NN_6 = cc(6,:)*tmp;
            
            A_tri_1 = ww'*NN_1';
            A_tri_2 = ww'*NN_2';
            A_tri_3 = ww'*NN_3';
            A_tri_4 = ww'*NN_4';
            A_tri_5 = ww'*NN_5';
            A_tri_6 = ww'*NN_6';
            
            A_tri = [A_tri_1; A_tri_2; A_tri_3; A_tri_4; A_tri_5; A_tri_6];
            
            Area_n(tri_ii) = Area_n(tri_ii) + A_tri;
            
        end
        
    end
end
