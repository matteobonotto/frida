

% % tic
% % L_VI = fun_L_VI_mex(tri, ...
% %     nn,...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     n_G_source, ...
% %     P_G_source,...
% %     w_G_source,...
% %     n_G_target, ....
% %     P_G_target,...
% %     w_G_target,...
% %     shape_f, ...
% %     1);
% % toc


clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(10,[Inf,2])'
spec{7} = 'coder.typeof(10,[Inf,1])'
spec{8} = 'coder.typeof(1,[1,1])'
spec{9} = 'coder.typeof(10,[Inf,2])'
spec{10} = 'coder.typeof(10,[Inf,1])'
spec{11} = 'coder.typeof(10,[Inf,100],[0 1])'
spec{12} = 'coder.typeof(1,[1,1])'

command = '';

for ii=1:12
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_L_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)


% % fun_L_VI_mex_1.mexw64 % parfor ciclo più esterno
% % fun_L_VI_mex_2.mexw64 % parfor secondo ciclo



%%


% % [M_VI] = fun_M_act_VI(tri, ...
% %     nn,...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     n_G_target, ....
% %     P_G_target,...
% %     w_G_target,...
% %     shape_f, ...
% %     n_act, ...
% %     tri_act, ...
% %     nodes_act, ...
% %     keyreg_act, ...
% %     OPT_PARALLEL)


clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(10,[Inf,2])'
spec{7} = 'coder.typeof(10,[Inf,1])'
spec{8} = 'coder.typeof(10,[Inf,100],[0 1])'
spec{9} = 'coder.typeof(1,[1,1])'
spec{10} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{11} = 'coder.typeof(10,[Inf,2])'
spec{12} = 'coder.typeof(10,[Inf,1])'
spec{13} = 'coder.typeof(1,[1,1])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_M_act_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)




%%

% % tic
% % R_VI = fun_R_VI(tri, ...
% %     nn,...
% %     N_order, ...
% %     n_G_source, ...
% %     P_G_source,...
% %     w_G_source,...
% %     shape_f, ...
% %     eta);
% % toc

clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(1,[1,1])'
spec{5} = 'coder.typeof(10,[Inf,2])'
spec{6} = 'coder.typeof(10,[Inf,1])'
spec{7} = 'coder.typeof(10,[Inf,100],[0 1])'
spec{8} = 'coder.typeof(10,[Inf,1])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_R_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)



%%


% % [M_VI] = fun_M_act_VI(tri, ...
% %     nn,...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     n_G_target, ....
% %     P_G_target,...
% %     w_G_target,...
% %     shape_f, ...
% %     n_act, ...
% %     tri_act, ...
% %     nodes_act, ...
% %     keyreg_act, ...
% %     OPT_PARALLEL)


clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(10,[Inf,2])'
spec{7} = 'coder.typeof(10,[Inf,1])'
spec{8} = 'coder.typeof(10,[Inf,100],[0 1])'
spec{9} = 'coder.typeof(1,[1,1])'
spec{10} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{11} = 'coder.typeof(10,[Inf,2])'
spec{12} = 'coder.typeof(10,[Inf,1])'
spec{13} = 'coder.typeof(1,[1,1])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_M_act_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)


%%

% % G_flux_VI = fun_G_flux_VI(tri, ...
% %     nn,...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     n_G_source, ...
% %     P_G_source,...
% %     w_G_source,...
% %     PSISENS_R, ...
% %     PSISENS_Z, ...
% %     shape_f, ...
% %     1);

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])';
spec{2} = 'coder.typeof(1,[1,1])';
spec{3} = 'coder.typeof(1,[1,1])';
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])';
spec{5} = 'coder.typeof(1,[1,1])';
spec{6} = 'coder.typeof(10,[Inf,2])';
spec{7} = 'coder.typeof(10,[Inf,1])';
spec{8} = 'coder.typeof(10,[Inf,1])';
spec{9} = 'coder.typeof(10,[Inf,1])';
spec{10} = 'coder.typeof(10,[Inf,100],[0 1])';
spec{11} = 'coder.typeof(1,[1,1])';

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_G_flux_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)


%%

% % G_flux_VI = fun_G_Aphi_VI(tri, ...
% %     nn,...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     n_G_source, ...
% %     P_G_source,...
% %     w_G_source,...
% %     PSISENS_R, ...
% %     PSISENS_Z, ...
% %     shape_f, ...
% %     1);

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])';
spec{2} = 'coder.typeof(1,[1,1])';
spec{3} = 'coder.typeof(1,[1,1])';
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])';
spec{5} = 'coder.typeof(1,[1,1])';
spec{6} = 'coder.typeof(10,[Inf,2])';
spec{7} = 'coder.typeof(10,[Inf,1])';
spec{8} = 'coder.typeof(10,[Inf,1])';
spec{9} = 'coder.typeof(10,[Inf,1])';
spec{10} = 'coder.typeof(10,[Inf,100],[0 1])';
spec{11} = 'coder.typeof(1,[1,1])';

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_G_Aphi_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)


%%

% % [G_Br_VI,G_Bz_VI] = fun_G_B_VI(tri, ...
% %     nn,...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     n_G_source, ...
% %     P_G_source,...
% %     w_G_source,...
% %     BSENS_R, ...
% %     BSENS_Z, ...
% %     shape_f, ...
% %     1);

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])';
spec{2} = 'coder.typeof(1,[1,1])';
spec{3} = 'coder.typeof(1,[1,1])';
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])';
spec{5} = 'coder.typeof(1,[1,1])';
spec{6} = 'coder.typeof(10,[Inf,2])';
spec{7} = 'coder.typeof(10,[Inf,1])';
spec{8} = 'coder.typeof(10,[Inf,1])';
spec{9} = 'coder.typeof(10,[Inf,1])';
spec{10} = 'coder.typeof(10,[Inf,100],[0 1])';
spec{11} = 'coder.typeof(1,[1,1])';

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_G_B_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)




%%

tic
L_VI = fun_L_VI_stable_mex_c(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    degree_G_target, ...
    n_G_target, ....
    P_G_target);
toc
tic
L_VI = fun_L_VI_stable_mex_cpp(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    degree_G_target, ...
    n_G_target, ....
    P_G_target);
toc
tic
L_VI = fun_L_VI_stable_mex_mingw(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    degree_G_target, ...
    n_G_target, ....
    P_G_target);
toc


clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(10,[Inf,2])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(1,[1,1])'
spec{7} = 'coder.typeof(10,[Inf,2])'
spec{8} = 'coder.typeof(1,[1,1])'
spec{9} = 'coder.typeof(1,[1,1])'
spec{10} = 'coder.typeof(10,[Inf,2])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

mexcfg.TargetLang = 'C';

command_all = ['codegen fun_L_VI_stable -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)



%%


% % tic
% % [M_VI] = fun_M_act_VI_stable(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     degree_G_source_ext, ...
% %     degree_G_target, ...
% %     n_G_target, ....
% %     P_G_target,...
% %     n_act, ...
% %     tri_act, ...
% %     nodes_act, ...
% %     keyreg_act);
% % toc


clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(10,[Inf,2])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(1,[1,1])'
spec{7} = 'coder.typeof(1,[1,1])'
spec{8} = 'coder.typeof(10,[Inf,2])'
spec{9} = 'coder.typeof(1,[1,1])'
spec{10} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{11} = 'coder.typeof(10,[Inf,2])'
spec{12} = 'coder.typeof(10,[Inf,1])'

command = '';

for ii=1:size(spec,2)
% % for ii=1:9
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex');
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;
mexcfg

command_all = ['codegen fun_M_act_VI_stable -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)


%%

tic
R_VI = fun_R_VI_stable(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    shape_f, ...
    eta);
toc



clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(10,[Inf,2])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(1,[1,1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(10,[Inf,2])'
spec{7} = 'coder.typeof(10,[Inf,100],[0 1])'
spec{8} = 'coder.typeof(10,[Inf,1])'

command = '';

for ii=1:size(spec,2)
% % for ii=1:9
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex');
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;
mexcfg

command_all = ['codegen fun_R_VI_stable -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)

%%

% % G_flux_VI = fun_G_flux_VI_stable(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source,...
% %     PSISENS_R, ...
% %     PSISENS_Z);

clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(10,[Inf,2])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(1,[1,1])'
spec{7} = 'coder.typeof(10,[Inf,2])'
spec{8} = 'coder.typeof(1,[Inf,1])'
spec{9} = 'coder.typeof(1,[Inf,1])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_G_flux_VI_stable -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)

%%
% % tic
% % G_Aphi_VI = fun_G_Aphi_VI_stable(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     MatInd_nodes_tri,...
% %     degree_G_source, ...
% %     n_G_source, ...
% %     P_G_source,...
% %     PSISENS_R, ...
% %     PSISENS_Z);
% % toc

clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(10,[Inf,2])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(1,[1,1])'
spec{7} = 'coder.typeof(10,[Inf,2])'
spec{8} = 'coder.typeof(1,[Inf,1])'
spec{9} = 'coder.typeof(1,[Inf,1])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_G_Aphi_VI_stable -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)

%%
tic
[G_Br_VI,G_Bz_VI] = fun_G_B_VI_stable(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    BSENS_R, ...
    BSENS_Z);

toc

clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(10,[Inf,2])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,20],[0 1])'
spec{5} = 'coder.typeof(1,[1,1])'
spec{6} = 'coder.typeof(1,[1,1])'
spec{7} = 'coder.typeof(10,[Inf,2])'
spec{8} = 'coder.typeof(1,[Inf,1])'
spec{9} = 'coder.typeof(1,[Inf,1])'

command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_G_B_VI_stable -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)






%%
%%
%%

L_VI = fun_assemby_L_VI(tri, ...
    nn, ...
    n_G_source, ...
    n_G_target, ....
    P_G_target, ...
    Green_Mat_Gauss_Aphi, ...
    det_Jac, ...
    w_G_soruce_norm, ...
    W_r_source, ...
    w_G_target_norm, ...
    W_r_target);



clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(1,[1,1])'
spec{5} = 'coder.typeof(10,[Inf,2])'
spec{6} = 'coder.typeof(10,[Inf,Inf],[0 1])'
spec{7} = 'coder.typeof(1,[Inf,1])'
spec{8} = 'coder.typeof(10,[100,1],[1 0])'
spec{9} = 'coder.typeof(10,[100,100],[1 1])'
spec{10} = 'coder.typeof(10,[100,1],[1 0])'
spec{11} = 'coder.typeof(10,[100,100],[1 1])'


command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_assemby_L_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)




%%

M_VI = fun_assembly_M_VI(M_VI, ...
    ii, ...
    tri, ...
    n_G_target, ...
    P_G_target, ...
    w_G_target_norm, ...
    vec_Aphi_all, ...
    N_nodes_sf, ...
    W_r_target, ...
    det_Jac);



clear spec

spec{1} = 'coder.typeof(10,[Inf,Inf],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{4} = 'coder.typeof(1,[1,1])'
spec{5} = 'coder.typeof(10,[Inf,2])'
spec{6} = 'coder.typeof(1,[Inf,1])'
spec{7} = 'coder.typeof(1,[Inf,1])'
spec{8} = 'coder.typeof(1,[1,1])'
spec{9} = 'coder.typeof(10,[100,100],[1 1])'
spec{10} = 'coder.typeof(1,[Inf,1])'



command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_assembly_M_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)



%%


G_flux_VI = fun_assemby_G_VI(tri, ...
    nn, ...
    n_G_source, ...
    Green_Mat_Gauss_flux, ...
    det_Jac, ...
    w_G_soruce_norm, ...
    W_r_source);



clear spec

spec{1} = 'coder.typeof(10,[Inf,10],[0 1])'
spec{2} = 'coder.typeof(1,[1,1])'
spec{3} = 'coder.typeof(1,[1,1])'
spec{4} = 'coder.typeof(10,[Inf,Inf],[0 1])'
spec{5} = 'coder.typeof(1,[Inf,1])'
spec{6} = 'coder.typeof(1,[Inf,1])'
spec{7} = 'coder.typeof(10,[100,100],[1 1])'


command = '';

for ii=1:size(spec,2)
    
    command = [command num2str(spec{ii}) ','];
    
end

clc
mexcfg = coder.config('mex')
mexcfg.IntegrityChecks = false;
mexcfg.ExtrinsicCalls = false;
mexcfg.ResponsivenessChecks = false;

command_all = ['codegen fun_assemby_G_VI -args {' command(1:end-1) '} -O enable:inline -O enable:openmp']

eval(command_all)






