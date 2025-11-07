perturbation = SETTINGS.perturbation;

% % fprintf('Computing Jacobian matrix (semi-analytic) \n')
calc_Jacobian_FRIDA_semianalytic_fast

DRpsibc_Dpsi = S_tilde;
DRpsibc_Dpsibc = S_tilde_bc;

if SETTINGS.FAR_FIELD_SPARS == true
    threshold = SETTINGS.FAR_FIELD_SPARS_THRESHOLD;
    DRpsibc_Dpsi_spars = fun_far_field_sparsification(DRpsibc_Dpsi,threshold);
else
    DRpsibc_Dpsi_spars = DRpsibc_Dpsi;
end

DR_psi   = [DRpsi_Dpsi DRpsi_Dpsibc DRpsi_Dpsia DRpsi_Dpsib DRpsi_Dlambda];
DRpsibc  = [DRpsibc_Dpsi_spars DRpsibc_Dpsibc DRpsibc_Dpsia DRpsibc_Dpsib DRpsibc_Dlambda];
DRpsia   = [sigma_a DRpsia_Dpsia DRpsia_Dpsib DRpsia_Dlambda];
DRpsib   = [sigma_b DRpsib_Dpsia DRpsib_Dpsib DRpsib_Dlambda];
DRlambda = [DRlambda_Dpsi DRlambda_Dpsibc DRlambda_Dpsia DRlambda_Dpsib DRlambda_Dlambda];

JAC = [DR_psi; ...
    DRpsibc; ...
    DRpsia; ...
    DRpsib; ...
    DRlambda];