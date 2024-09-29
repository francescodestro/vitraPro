##############################################################################
#
#   Function for defining the model parameters 
#
###############################################################################
#
#   Detailed description of the parameters: Table S2 in Destro and Braatz, 2024 
#   For STV/DIP systems: STV = virus 1, DIP = virus 2.
#   For systems with two STVs, no predefined order.
#
###############################################################################

import numpy as np

def parameters(age_viab, age_end_inf, Dtau, scaling):

    # Cell growth and death
    mu_T = 0.028  # 1/h
    mu_I1 = 0  # 1/h
    k_deathT = 8e-5  # 1/h
    tau_death_I1 = 24  # hpi
    k_death_I1 = 0.0029  # 1/h

    # Viral binding
    k_bind0_I1 = 8e-7 * scaling  # mL/cell/h
    tau0_bind_I1 = 1.80  # hpi
    beta_bind_I1 = 0.5  # –

    # Virus internalization
    k_int_I1 = 0.01 * 60  # 1/h
    eta_I1 = 0.5  # -

    # Viral replication
    k_repl_I1 = 0.7318  # 1/h
    replStart_I1 = 6  # hpi
    replEnd_I1 = 18  # hpi

    # Viral progeny release
    k_rel_I1 = 5  # PFU/cell/h
    tau_rel_on_I1 = 18  # hpi
    tau_rel_off_I1 = 72  # hpi

    # Viral genome degradation in nucleus
    k_d_N1 = 0  # 1/h

    # Viral degradation in supernatant
    k_d_V1 = 7e-3  # 1/h

    # Substrate
    K_s = 1.3 * 1e3  # nmol/mL
    Ys_T = 1.2e-4  # nmol/cell/h
    Ys_I1 = 0  # nmol/cell/h

    # Nonviable cell lysis
    k_lys = 0  # 1/h


    ############################################################################
    #
    # don't modify part below
    #
    ############################################################################
    
    # Mesh and state vectors indexing
    n_bins = int(age_viab / Dtau)
    n_bins_inf = int(age_end_inf / Dtau)
    age_bins = np.arange(0, age_viab, Dtau)

    # First state vector
    startI1 = 2
    endI1 = 1 + n_bins

    # Other state vectors: internalized virus
    startB1 = 0
    endB1 = n_bins - 1

    # Adjust substrate concentration
    K_s = K_s / scaling

    # Parameters for infected cells
    f_rel_I1 = np.zeros(len(age_bins))
    f_rel_I1[age_bins >= tau_rel_on_I1] = 1
    f_rel_I1[age_bins >= tau_rel_off_I1] = 0
    f_rel_I1[-1] = 0

    f_repl_I1 = (0.5 + 0.5 * np.tanh((age_bins - replStart_I1) / 0.3)) * (0.5 + 0.5 * np.tanh((replEnd_I1 - age_bins) / 1))

    n1 = 0.5 + 0.5 * np.tanh((age_bins[:n_bins_inf] - tau0_bind_I1) / 0.01)
    k_bind_I1 = k_bind0_I1 * ((1 - n1) + n1 * np.exp(-beta_bind_I1 * (age_bins[:n_bins_inf] - tau0_bind_I1)))

    # Correction to reallocation when age_end_inf == age_viab
    a = -1 if age_end_inf == age_viab else 0

    # Solver initialization parameters
    hmax0 = Dtau  # CFL condition
    rtol = 1.e-6
    atol = 1.e-6
    threshold = atol / rtol

    # Prepare parameters dictionary to pass to solver
    p = {
        "startI1": startI1,
        "endI1": endI1,
        "startB1": startB1,
        "endB1": endB1,
        "age_end_inf": age_end_inf,
        "a": a,
        "n_bins": n_bins,
        "n_bins_inf": n_bins_inf,
        "age_bins": age_bins,
        "tau_death_I1": tau_death_I1,
        "k_death_I1": k_death_I1,
        "k_bind_I1": k_bind_I1,
        "mu_T": mu_T,
        "mu_I1": mu_I1,
        "k_deathT": k_deathT,
        "k_int_I1": k_int_I1,
        "eta_I1": eta_I1,
        "f_repl_I1": f_repl_I1,
        "k_repl_I1": k_repl_I1,
        "f_rel_I1": f_rel_I1,
        "k_rel_I1": k_rel_I1,
        "k_degrBV_1": k_d_V1,
        "k_lys": k_lys,
        "k_d_N1": k_d_N1,
        "K_s": K_s,
        "Ys_T": Ys_T,
        "Ys_I1": Ys_I1,
        "threshold": threshold,
        "rtol": rtol,
        "Dtau": Dtau,
        "hmax0": hmax0
    }

    return p
