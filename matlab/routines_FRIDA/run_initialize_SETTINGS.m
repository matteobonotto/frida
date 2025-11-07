%%

% % warning('New release available')



%% Initialization of SETTINGS 

if ~any(strcmp('PREPROC',fieldnames(SETTINGS))); SETTINGS.PREPROC = true; end

if ~any(strcmp('POSTPROC',fieldnames(SETTINGS))); SETTINGS.POSTPROC = false; end

if ~any(strcmp('SAVE_OUTPUT',fieldnames(SETTINGS))); SETTINGS.SAVE_OUTPUT = false; end

if ~any(strcmp('PREPROC_REORDERING',fieldnames(SETTINGS))); SETTINGS.PREPROC_REORDERING = false; end

if ~any(strcmp('VAC_METHOD',fieldnames(SETTINGS))); SETTINGS.VAC_METHOD = 2; end

if ~any(strcmp('RUN_MEX_ROUTINE',fieldnames(SETTINGS))); SETTINGS.RUN_MEX_ROUTINE = false; end

if ~any(strcmp('FIGURES',fieldnames(SETTINGS))); SETTINGS.FIGURES = 1; end

if ~any(strcmp('FIGURES_DEBUG',fieldnames(SETTINGS))); SETTINGS.FIGURES_DEBUG = 0; end

if ~any(strcmp('GAUSS_QUAD_DEGREE_PLA',fieldnames(SETTINGS))); SETTINGS.GAUSS_QUAD_DEGREE_PLA = 10; end % 24 Gauss Points

if ~any(strcmp('GAUSS_QUAD_DEGREE_STIFFMAT',fieldnames(SETTINGS))); SETTINGS.GAUSS_QUAD_DEGREE_STIFFMAT = 8; end % 16 Gauss Points

if ~any(strcmp('N_POINTS_SEPARATRIX',fieldnames(SETTINGS))); SETTINGS.N_POINTS_SEPARATRIX = 200; end % Points for separatrix computation

if ~any(strcmp('TOLL',fieldnames(SETTINGS))); SETTINGS.TOLL = 1e-7; end

if ~any(strcmp('SCARTO',fieldnames(SETTINGS))); SETTINGS.SCARTO = 1e-7; end

if ~any(strcmp('ii_freeB_max',fieldnames(SETTINGS))); SETTINGS.ii_freeB_max = 20; end

if ~any(strcmp('THRESHOLD_FIXB',fieldnames(SETTINGS))); SETTINGS.THRESHOLD_FIXB = 1E-3; end

if ~any(strcmp('TOLL_FIXB',fieldnames(SETTINGS))); SETTINGS.TOLL_FIXB = 5E-4; end

if ~any(strcmp('N_ITERMAX_FIXB',fieldnames(SETTINGS))); SETTINGS.N_ITERMAX_FIXB = 20; end

if ~any(strcmp('NFREEB_FIXB',fieldnames(SETTINGS))); SETTINGS.NFREEB_FIXB = 10; end

if ~any(strcmp('perturbation',fieldnames(SETTINGS))); SETTINGS.perturbation = sqrt(eps); end

if ~any(strcmp('SOLVER_RELAXED',fieldnames(SETTINGS))); SETTINGS.SOLVER_RELAXED = true; end

if ~any(strcmp('NN_SIZE_MAX_PARALLEL_JACOBIAN',fieldnames(SETTINGS))); SETTINGS.NN_SIZE_MAX_PARALLEL_JACOBIAN = 100; end

if ~any(strcmp('NK_preconditioner',fieldnames(SETTINGS))); SETTINGS.NK_preconditioner = 1; end

if ~any(strcmp('SMOOTH_FDFPROFILE',fieldnames(SETTINGS))); SETTINGS.SMOOTH_FDFPROFILE = false; end

if ~any(strcmp('SMOOTH_DPPROFILE',fieldnames(SETTINGS))); SETTINGS.SMOOTH_DPPROFILE = false; end

if ~any(strcmp('FAR_FIELD_SPARS',fieldnames(SETTINGS))); SETTINGS.FAR_FIELD_SPARS = true; end

if ~any(strcmp('FAR_FIELD_SPARS',fieldnames(SETTINGS))); SETTINGS.FAR_FIELD_SPARS = true; end

if ~any(strcmp('FAR_FIELD_SPARS_THRESHOLD',fieldnames(SETTINGS))); SETTINGS.FAR_FIELD_SPARS_THRESHOLD = .005; end

if ~isfield(SETTINGS,'QuasiNewton'); SETTINGS.QuasiNewton = false; end

if ~isfield(SETTINGS,'QuasiNewton_factor'); SETTINGS.QuasiNewton_factor = 0.2; end

if ~isfield(SETTINGS,'VERBOSE'); SETTINGS.VERBOSE = true; end

if ~isfield(SETTINGS,'IS_EVOL'); SETTINGS.IS_EVOL = true; end

if ~isfield(SETTINGS,'RUN_FIXBOUNDARY'); SETTINGS.RUN_FIXBOUNDARY = true; end

if ~isfield(SETTINGS,'INITIAL_PLASMA_MODEL'); SETTINGS.INITIAL_PLASMA_MODEL = 1; end

if ~isfield(SETTINGS,'ii_START'); SETTINGS.ii_START = 1; end



if ~any(strcmp('SOLVER_RELAX_TRESHOLD_1',fieldnames(SETTINGS)))
    
    SETTINGS.SOLVER_RELAX_TRESHOLD_1 = 30;
    SETTINGS.SOLVER_RELAX_TRESHOLD_2 = 50; 
    SETTINGS.SOLVER_RELAX_TRESHOLD_3 = 150;
    SETTINGS.SOLVER_RELAX_TRESHOLD_4 = 500;
    
    SETTINGS.SOLVER_RELAX_ALPHA_1 = 1;
    SETTINGS.SOLVER_RELAX_ALPHA_2 = .5;
    SETTINGS.SOLVER_RELAX_ALPHA_3 = .2;
    SETTINGS.SOLVER_RELAX_ALPHA_4 = .075;
    SETTINGS.SOLVER_RELAX_ALPHA_5 = .01;
    
end

if ~any(strcmp('SOLVER',fieldnames(SETTINGS)))

    SETTINGS.SOLVER = 'NR';
    
    if ~any(strcmp('FAR_FIELD_SPARS',fieldnames(SETTINGS))); SETTINGS.FAR_FIELD_SPARS = true; end
    if ~any(strcmp('FAR_FIELD_SPARS_THRESHOLD',fieldnames(SETTINGS))); SETTINGS.FAR_FIELD_SPARS_THRESHOLD = .3; end
        
end


if any(strcmp('alpha_M',fieldnames(plaparameter)))
    SETTINGS.J_PARAMETRIZATION_TYPE = 2;
else 
    SETTINGS.J_PARAMETRIZATION_TYPE = 1;
end


%% Suppress run of MEX files if not on WINDOWS environment

if ismac
    SETTINGS.RUN_MEX_ROUTINE = false;
    
elseif isunix
    SETTINGS.RUN_MEX_ROUTINE = false;
    
end


%% INPUT

if any(strcmp('plaparameter',fieldnames(SETTINGS)))
    
    % modify Centroid
    if any(strcmp('Centroid',fieldnames(SETTINGS.plaparameter)))
        plaparameter.Centroid = SETTINGS.plaparameter.Centroid;
    end
    
    % modify Ipla
    if any(strcmp('Ipla',fieldnames(SETTINGS.plaparameter)))
        plaparameter.Ipla = SETTINGS.plaparameter.Ipla;
    end
    
    % modify beta_0
    if any(strcmp('beta_0',fieldnames(SETTINGS.plaparameter)))
        plaparameter.beta_0 = SETTINGS.plaparameter.beta_0;
    end
    
    % modify R_0
    if any(strcmp('R_0',fieldnames(SETTINGS.plaparameter)))
        plaparameter.R_0 = SETTINGS.plaparameter.R_0;
    end
    
    % modify alpha_M
    if any(strcmp('alpha_M',fieldnames(SETTINGS.plaparameter)))
        plaparameter.alpha_M = SETTINGS.plaparameter.alpha_M;
    end
    
    % modify alpha_N
    if any(strcmp('alpha_N',fieldnames(SETTINGS.plaparameter)))
        plaparameter.alpha_N = SETTINGS.plaparameter.alpha_N;
    end
    
    % modify psibar
    if any(strcmp('psibar',fieldnames(SETTINGS.plaparameter)))
        plaparameter.psibar = SETTINGS.plaparameter.psibar;
    end
    
    % modify FdF
    if any(strcmp('FdF',fieldnames(SETTINGS.plaparameter)))
        plaparameter.FdF = SETTINGS.plaparameter.FdF;
    end
    
    % modify dP
    if any(strcmp('dP',fieldnames(SETTINGS.plaparameter)))
        plaparameter.dP = SETTINGS.plaparameter.dP;
    end
    
end















