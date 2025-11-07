function [RES] = fun_ResidueFRIDA_evol(INPUT_Residue)


%%
xx_k         = INPUT_Residue.xx_k;
meshData_loc = INPUT_Residue.meshData_loc;
solk         = INPUT_Residue.solk;
plaparameter = INPUT_Residue.plaparameter;
AA_rid       = INPUT_Residue.AA_rid;
Conductors   = INPUT_Residue.Conductors;
SETTINGS     = INPUT_Residue.SETTINGS;


%%
mu0=4*pi*1.e-7;

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    Centroid = plaparameter.Centroid;
    Ipla =     plaparameter.Ipla;
    psibar =   plaparameter.psibar;
    FdF =      plaparameter.FdF;
    dP =      plaparameter.dP;
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    R_0 =      plaparameter.R_0;
    beta_0 =   plaparameter.beta_0;
    alpha_M =  plaparameter.alpha_M;
    alpha_N =  plaparameter.alpha_N;
    Ipla =     plaparameter.Ipla;
end


nn=meshData_loc.nn;

nn_D = length(meshData_loc.ind_D);
nn_B = length(meshData_loc.ind_B);

ind_D = meshData_loc.ind_D;
ind_B = meshData_loc.ind_B;

rr_Gauss = meshData_loc.nodes_Gauss_pla(:,1);

%%
tri             = meshData_loc.t;
N_order         = meshData_loc.shape_functions_N_order;
P_Gauss         = meshData_loc.nodes_Gauss_pla;
n_Gauss         = meshData_loc.n_Gauss;
shape_functions = meshData_loc.shape_functions;

Psi_k=solk.Psi;
Psi_axis=solk.Psi_axis_node;
Psi_Bpla=sum(solk.Psi_Bpla)./numel(solk.Psi_Bpla);
psi_bar=(Psi_k-Psi_axis)/(Psi_Bpla-Psi_axis); %normalized poloidal magnetic flux

if SETTINGS.RUN_MEX_ROUTINE == false
    [psi_bar] = fun_evaluateFluxGaussPoints_v2(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
elseif SETTINGS.RUN_MEX_ROUTINE == true
    [psi_bar] = fun_evaluateFluxGaussPoints_v2_mex(tri,psi_bar,N_order,P_Gauss,n_Gauss,shape_functions);
end
solk.psi_bar = psi_bar;

if SETTINGS.J_PARAMETRIZATION_TYPE == 1
    
% %     psi_bar(psi_bar>1) = NaN;
    psi_bar(psi_bar<0) = 0;

    fdfn=interp1(psibar,FdF,psi_bar,'linear');   
    dpn=interp1(psibar,dP,psi_bar,'linear');   
    gg_tot=(rr_Gauss.*dpn+fdfn./(mu0*rr_Gauss)); %computation of nodes plasma current density
    gg_tot=gg_tot/2/pi;
    
    qq = isnan(gg_tot);
    gg_tot(isnan(gg_tot)) = 0;
    
    nn = meshData_loc.nn;
    ind_t_selected = meshData_loc.ind_t_INpla;
    ww_Gauss = meshData_loc.ww_Gauss_pla;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
        [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
    
    % %     figure;  plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),Source,'o')
    
elseif SETTINGS.J_PARAMETRIZATION_TYPE == 2
    
    tmp_gg = ((1-psi_bar.^alpha_M).^alpha_N);
    tmp_gg(imag(tmp_gg) ~= 0) = abs(tmp_gg(imag(tmp_gg) ~= 0)); % occhio!!!
    aa_r_tot = (rr_Gauss*beta_0/R_0+R_0*(1-beta_0)./rr_Gauss);
    gg_tot = tmp_gg.*aa_r_tot;
    
    nn = meshData_loc.nn;
    ind_t_selected = meshData_loc.ind_t_INpla;
    ww_Gauss = meshData_loc.ww_Gauss_pla;
    
    if SETTINGS.RUN_MEX_ROUTINE == false
        [Source] = fun_calcSourceFEM_v2(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    elseif SETTINGS.RUN_MEX_ROUTINE == true
        if size(ind_t_selected,1) < size(ind_t_selected,2); ind_t_selected = ind_t_selected'; end
        [Source] = fun_calcSourceFEM_v2_mex(tri,nn,gg_tot,ind_t_selected,N_order,P_Gauss,ww_Gauss,n_Gauss,shape_functions);
    end
    
end

Source_D = Source(ind_D);



%% Build stiffness matrix and RHS


AA_I = sparse(size(AA_rid,1),size(AA_rid,1));
AA_I(1:nn_D,nn+3) = -Source_D;
AA_I(nn+3,nn+3) = sum(Source_D);

AA = AA_rid + AA_I;

%%

RES = AA*xx_k - solk.RHS;

% % qq = load('temp_NK.mat')
% % 
% % norm(RES - qq.RR_psi_k)
% % norm(xx_k - qq.xx_k)
% % norm(AA*xx_k - qq.AA*qq.xx_k)
% % norm(RES - qq.RR_psi_k)



% % qq = [AA*xx_k RHS]



% % xx_k1 = AA\RHS;
% % 
% % Psi_tot = zeros(meshData_loc.nn,1);
% % Psi_tot(ind_D) = xx_k1(1:nn_D);
% % Psi_tot(ind_B) = xx_k1(nn_D+1:nn);
% % 
% % % % figure
% % % % plot3(meshData_loc.n(:,1),meshData_loc.n(:,2),Psi_tot,'.')
% % 
% % 
% % lambda_norm = xx_k1(end);
% % lambda = lambda_norm/mu0;
% % Psi_a = xx_k1(end-2);
% % Psi_b = xx_k1(end-1)
% % 
% % 
% % figure
% % axis equal;
% % Rplot=[.9*min(meshData_loc.n(:,1)) 1.08*max(meshData_loc.n(:,1))];
% % Zplot=1.3*[min(meshData_loc.n(:,2)) max(meshData_loc.n(:,2))];
% % rgrid=linspace(Rplot(1),Rplot(2),200);
% % zgrid=linspace(Zplot(1),Zplot(2),200);
% % [RR,ZZ] = meshgrid(rgrid,zgrid);
% % axis([Rplot Zplot]),  hold on,
% % PSI = griddata(meshData_loc.n(:,1),meshData_loc.n(:,2),Psi_tot,RR,ZZ);
% % contour(RR,ZZ,PSI,50);
% % colormap('cool'); colorbar vert;
% % contour(RR,ZZ,PSI,[Psi_b Psi_b],'k');
% % plot(meshData_loc.n(solk.ind_n_axis,1),meshData_loc.n(solk.ind_n_axis,2),'r*','LineWidth',2)
% % plot(meshData_loc.n(solk.ind_n_XP,1),meshData_loc.n(solk.ind_n_XP,2),'b*','LineWidth',2)
% % plot(meshData_loc.n(meshData_loc.ind_n_FW,1),meshData_loc.n(meshData_loc.ind_n_FW,2),'k-')







