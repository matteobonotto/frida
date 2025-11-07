

%%

L_VI              = meshData_ext.L_VI;
R_VI              = meshData_ext.R_VI;
U_VI              = meshData_ext.U_VI;
V_VI              = meshData_ext.V_VI;
D_VI              = meshData_ext.D_VI;
M_VI              = meshData_ext.M_VI;
G_BCs_VI_act      = meshData_ext.G_BCs_VI_act;
G_BCs_VI_pas      = meshData_ext.G_BCs_VI_pas;
G_CS_VI_act       = meshData_ext.G_CS_VI_act;
G_CS_VI_pas       = meshData_ext.G_CS_VI_pas;
G_flux_VI_pas     = meshData_ext.G_flux_VI_pas;
G_B_pickup_VI_pas = meshData_ext.G_B_pickup_VI_pas;
G_flux_VI_act     = meshData_ext.G_flux_VI_act;
G_B_pickup_VI_act = meshData_ext.G_B_pickup_VI_act;
vec_A             = meshData_ext.vec_A;

T_Ipla_eq =meshData_pla.T_Ipla_eq;
T_Ipla_eq_D = T_Ipla_eq(:,meshData_pla.ind_D);
meshData_pla.T_Ipla_eq_D = T_Ipla_eq_D;


%%
h_step = Conductors.time_sim(2:end) - Conductors.time_sim(1:end-1);

if sum(h_step)/numel(h_step) - h_step(1) > 1e-7
    warning('WARNING: non constant time steps! Choosing h_step = h_step(1).')
end
h_step = h_step(1);

M_act_pas = M_VI;
M_CS_pas = G_CS_VI_pas;
% % meshData_ext.M_CS_pas = M_CS_pas;

ind_act = (1:size(M_act_pas,2))';
ind_CS = (ind_act(end)+1:size(M_act_pas,2)+size(M_CS_pas,2))';

D_inc = zeros(size(D_VI));
D_inc(find(D_VI)) = 1;

ind_passive_cut = meshData_ext.ind_passive_cut;
L_VI_inv = L_VI\eye(size(L_VI));
beta_A = vec_A(:,ind_passive_cut).'*L_VI_inv;

H_act = beta_A*M_act_pas;
H_pla = beta_A*M_CS_pas*T_Ipla_eq_D;

A_SS_ODE = -R_VI*L_VI_inv;

EE = blkdiag(eye(size(L_VI)),zeros(numel(ind_passive_cut)));
A = [A_SS_ODE  V_VI*D_inc(:,ind_passive_cut); ...
    beta_A zeros(length(ind_passive_cut))];

B_act = [-A_SS_ODE*M_act_pas; -H_act];
B_pla = [-A_SS_ODE*M_CS_pas*T_Ipla_eq_D; -H_pla];

C = [L_VI_inv zeros(size(L_VI_inv,1),length(ind_passive_cut))];

D_act = -L_VI_inv*M_act_pas;
D_pla = -L_VI_inv*M_CS_pas*T_Ipla_eq_D;

meshData_ext.E = EE;
meshData_ext.A = A;
meshData_ext.C = C;
meshData_ext.B_act = B_act;
meshData_ext.B_pla = B_pla;
meshData_ext.D_act = D_act;
meshData_ext.D_pla = D_pla;
% % 
% % meshData_ext.L_ode_inv = L_VI_inv;


%%% For Eddy current equation
KK_D  = meshData_pla.KK_nobc(meshData_pla.ind_D,meshData_pla.ind_D);
KK_bc = meshData_pla.KK_nobc(meshData_pla.ind_D,meshData_pla.ind_B);

A_CN_1 = (EE - h_step*A/2);
A_CN_1_inv = A_CN_1\eye(size(A_CN_1));

A_CN_2 = (EE + h_step*A/2);

Q_tilde = h_step*.5*A_CN_1_inv*B_pla*KK_D/mu0;
Q_tilde_bc = h_step*.5*A_CN_1_inv*B_pla*KK_bc/mu0;

GG_tilde_bc_pas = G_BCs_VI_pas/2/pi;

P_tilde = -GG_tilde_bc_pas*( C*Q_tilde + D_pla*KK_D/mu0);
P_tilde_bc = -GG_tilde_bc_pas*( C*Q_tilde_bc + D_pla*KK_bc/mu0);

meshData_ext.h_step = h_step;
meshData_ext.A_CN_1 = A_CN_1;
meshData_ext.A_CN_1_inv = A_CN_1_inv;
meshData_ext.A_CN_2 = A_CN_2;
meshData_ext.P_tilde = P_tilde;
meshData_ext.P_tilde_bc = P_tilde_bc;
meshData_ext.Q_tilde = Q_tilde;
meshData_ext.Q_tilde_bc = Q_tilde_bc;
meshData_ext.GG_tilde_bc_pas = GG_tilde_bc_pas;



