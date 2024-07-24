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

def parameters(age_viab, age_end_inf, Dtau, scaling, DIP):
    
    # Cell growth and death
    mu_T = 0.028  # 1/h
    mu_I1 = 0  # 1/h
    mu_I2 = 0  # 1/h
    mu_Co = 0  # 1/h
    k_deathT = 8e-5  # 1/h
    tau_death_I1 = 24  # hpi
    tau_death_I2 = tau_death_I1
    tau_death_Co = tau_death_I1
    k_death_I1 = 0.0029  # 1/h
    k_death_I2 = k_death_I1  # 1/h
    k_death_Co = k_death_I1  # 1/h
    f_death_Co1 = 1  # 1/h
    f_death_Co2 = f_death_Co1  # 1/h

    # Viral binding
    k_bind0_I1 = 8e-7 * scaling  # mL/cell/h
    tau0_bind_I1 = 1.80  # hpi
    beta_bind_I1 = 0.5  # -
    
    k_bind0_I2 = k_bind0_I1  # mL/cell/h
    tau0_bind_I2 = tau0_bind_I1  # hpi
    beta_bind_I2 = beta_bind_I1  # -
    
    tau0_bind_Co1 = tau0_bind_I1  # hpi
    beta_bind_Co1 = beta_bind_I1  # -
    
    tau0_bind_Co2 = tau0_bind_I1  # hpi
    beta_bind_Co2 = beta_bind_I1  # -
    
    # Virus internalization
    k_int_I1 = 0.01 * 60  # 1/h
    k_int_I2 = k_int_I1  # 1/h
    k_int_Co1 = k_int_I1  # 1/h
    k_int_Co2 = k_int_I1  # 1/h
    eta_I1 = 0.5  # -
    eta_I2 = eta_I1  # -
    eta_Co1 = eta_I1  # -
    eta_Co2 = eta_I1  # -
    
    # Viral replication
    k_repl_I1 = 0.7318  # 1/h
    replStart_I1 = 6  # hpi
    replEnd_I1 = 18  # hpi
    k_repl_I2 = k_repl_I1  # 1/h
    replStart_I2 = replStart_I1  # hpi
    replEnd_I2 = replEnd_I1  # hpi
    k_repl_Co1 = k_repl_I1  # 1/h
    replStart_Co1 = replStart_I1  # hpi
    replEnd_Co1 = replEnd_I1  # hpi
    k_repl_Co2 = k_repl_I1  # 1/h
    replStart_Co2 = replStart_I1  # hpi
    replEnd_Co2 = replEnd_I1  # hpi
    
    # Viral progeny release
    k_rel_I1 = 5  # PFU/cell/h
    tau_rel_on_I1 = 18  # hpi
    tau_rel_off_I1 = 72  # hpi
    k_rel_I2 = k_rel_I1  # PFU/cell/h
    tau_rel_on_I2 = tau_rel_on_I1  # hpi
    tau_rel_off_I2 = tau_rel_off_I1  # hpi
    k_rel_Co = k_rel_I1  # PFU/cell/h
    tau_rel_on_Co = tau_rel_on_I1  # hpi
    tau_rel_off_Co = tau_rel_off_I1  # hpi
    randV2_I1 = 0  # -
    randV1_I2 = 0  # -
    
    # Viral genome degradation in nucleus
    k_d_N1 = 0  # 1/h
    k_d_N2 = 0  # 1/h
    
    # Viral degradation in supernatant
    k_d_V1 = 7e-3  # 1/h
    k_d_V2 = k_d_V1  # 1/h
    
    # Substrate
    K_s = 1.3 * 1e3  # nmol/mL
    Ys_T = 1.2e-4  # nmol/cell/h
    Ys_I1 = 0  # nmol/cell/h
    Ys_I2 = 0  # nmol/cell/h
    Ys_C = 0  # nmol/cell/h
    
    # Nonviable cell lysis
    k_lys=0 # 1/h

    #########################################################################
    # Don't edit part below!
    #########################################################################

    # Mesh and state vectors indexing
    # Mesh
    n_bins = int(age_viab / Dtau)
    n_bins_inf = int(age_end_inf / Dtau)
    age_bins = np.arange(0, age_viab, Dtau)

    # First state vector
    startI1 = 3
    endI1 = 2 + n_bins
    startI2 = 3 + n_bins
    endI2 = 2 + 2 * n_bins
    startCo = 3 + 2 * n_bins
    endCo = 2 + 2 * n_bins + n_bins ** 2 - (n_bins - n_bins_inf) * (n_bins - n_bins_inf + 1)
    
    # Other state vectors: internalized virus
    startB1 = 0
    endB1 = n_bins-1
    startB2 = n_bins
    endB2 = 2 * n_bins -1
    startB1Co = 2 * n_bins
    endB1Co = 2 * n_bins + n_bins ** 2 - (n_bins - n_bins_inf) * (n_bins - n_bins_inf + 1)-1
    startB2Co = endB1Co + 1
    endB2Co = endB1Co + n_bins ** 2 - (n_bins - n_bins_inf) * (n_bins - n_bins_inf + 1)
    
    ind = np.zeros(n_bins ** 2,dtype=int)
    ind_scCo = np.zeros(n_bins ** 2,dtype=int)
    ind_scB1Co = np.zeros(n_bins ** 2,dtype=int)
    ind_scB2Co = np.zeros(n_bins ** 2,dtype=int)
    
    K_s /= scaling

    # Parameters for coinfected cells
    k_bind0_Co1 = k_bind0_I1  # mL/cell/h
    k_bind0_Co2 = k_bind0_I1  # mL/cell/h
    
    k_bind_Co1 = np.zeros(n_bins ** 2)
    k_bind_Co2 = np.zeros(n_bins ** 2)
    
    age_co = np.zeros(n_bins ** 2)
    f_repl_Co1 = np.zeros(n_bins ** 2)
    f_repl_Co2 = np.zeros(n_bins ** 2)

    i = 0
    for i1 in range(0, n_bins):
        for i2 in range(0, n_bins):
            if (i1 - i2) ** 2 < n_bins_inf ** 2:
                ind[i2 + n_bins * i1] = int(i)+1  # entry of dCodt, age_co, f_repl, f_rel
                ind_scCo[i2 + n_bins * (i1)] =  int(i) + startCo   # entry of dxdt
                ind_scB1Co[i2 + n_bins * (i1)] =  int(i) + startB1Co   # entry of dbind_dt and dnucl_dt
                ind_scB2Co[i2 + n_bins * (i1)] =  int(i) + startB2Co   # entry of dbind_dt and dnucl_dt
                
                if DIP:  # infection age for STV/DIP co-infected cells
                    current_age = (i1) * Dtau
                    if i2 == n_bins:
                        current_age = (i2) * Dtau
                else:  # infection age for STV1/STV2 co-infected cells
                    current_age = (max(i1, i2) ) * Dtau
                
                n1 = 0.5 + 0.5 * np.tanh((current_age - tau0_bind_Co1) / 0.01)
                k_bind_Co1[i2 + n_bins * (i1 )] = k_bind0_Co1 * ((1 - n1) + n1 * np.exp(-beta_bind_Co1 * (current_age - tau0_bind_Co1)))
                
                n2 = 0.5 + 0.5 * np.tanh((current_age - tau0_bind_Co2) / 0.01)
                k_bind_Co2[i2 + n_bins * (i1 )] = k_bind0_Co2 * ((1 - n1) + n2 * np.exp(-beta_bind_Co2 * (current_age - tau0_bind_Co2)))
                
                f_repl_Co1[i2 + n_bins * (i1 )] = (0.5 + 0.5 * np.tanh((current_age - replStart_Co1) / 0.3)) * \
                    (0.5 + 0.5 * np.tanh((replEnd_Co1 - current_age) / 1))
                f_repl_Co2[i2 + n_bins * (i1 )] = (0.5 + 0.5 * np.tanh((current_age - replStart_Co2) / 0.3)) * \
                    (0.5 + 0.5 * np.tanh((replEnd_Co2 - current_age) / 1))
                
                age_co[i2 + n_bins * (i1)] = current_age
                i += 1
                
    k_bind_Co1[age_co>age_end_inf]=0
    k_bind_Co2[age_co>age_end_inf]=0
            
    age_co = age_co[ind > 0]
    f_repl_Co1 = f_repl_Co1[ind > 0]
    f_repl_Co2 = f_repl_Co2[ind > 0]
    k_bind_Co1 = k_bind_Co1[ind > 0]
    k_bind_Co2 = k_bind_Co2[ind > 0]
    
    ind=ind-1

    f_rel_Co = np.zeros(len(age_co))
    f_rel_Co[age_co >= tau_rel_on_Co] = 1
    f_rel_Co[age_co <= tau_rel_on_Co] = 0
    f_rel_Co[age_co >= tau_rel_off_Co] = 0
    f_rel_Co[age_co == (age_viab - Dtau)] = 0

    f_rel_I1 = np.zeros(len(age_bins))
    f_rel_I1[age_bins >= tau_rel_on_I1] = 1
    f_rel_I1[age_bins >= tau_rel_off_I1] = 0
    f_rel_I1[-1] = 0
    
    f_rel_I2 = np.zeros(len(age_bins))
    f_rel_I2[age_bins >= tau_rel_on_I2] = 1
    f_rel_I2[age_bins >= tau_rel_off_I2] = 0
    f_rel_I2[-1] = 0
    
    f_repl_I1 = (0.5 + 0.5 * np.tanh((age_bins - replStart_I1) / 0.3)) * \
                (0.5 + 0.5 * np.tanh((replEnd_I1 - age_bins) / 1))
    
    f_repl_I2 = (0.5 + 0.5 * np.tanh((age_bins - replStart_I2) / 0.3)) * \
                (0.5 + 0.5 * np.tanh((replEnd_I2 - age_bins) / 1))
    
    n1 = 0.5 + 0.5 * (np.tanh((age_bins[:n_bins_inf] - tau0_bind_I1) / 0.01))
    k_bind_I1 = k_bind0_I1 * ((1 - n1) + n1 * np.exp(-beta_bind_I1 * (age_bins[:n_bins_inf] - tau0_bind_I1)))
    
    n2 = 0.5 + 0.5 * (np.tanh((age_bins[:n_bins_inf] - tau0_bind_I2) / 0.01))
    k_bind_I2 = k_bind0_I2 * ((1 - n2) + n2 * np.exp(-beta_bind_I2 * (age_bins[:n_bins_inf] - tau0_bind_I2)))
    
    # correction to reallocation when age_end_inf = age_viability
    if age_end_inf == age_viab:
        a = -1 
    else:
        a = 0
    
    # Initialization of solver and state vector
    hmax0 = Dtau  # CFL condition
    rtol = 1.e-6
    atol = 1.e-6
    threshold = atol / rtol
    
    p = {
        "startI1": startI1,
        "endI1": endI1,
        "startI2": startI2,
        "endI2": endI2,
        "startCo": startCo,
        "endCo": endCo,
        "startB1": startB1,
        "startB2": startB2,
        "endB1": endB1,
        "endB2": endB2,
        "endB1Co": endB1Co,
        "endB2Co": endB2Co,
        "startB1Co": startB1Co,
        "startB2Co": startB2Co,
        "ind_scCo": ind_scCo,
        "age_end_inf": age_end_inf,
        "a": a,
        "ind": ind,
        "ind_scB1Co": ind_scB1Co,
        "ind_scB2Co": ind_scB2Co,
        "n_bins": n_bins,
        "n_bins_inf": n_bins_inf,
        "age_bins": age_bins,
        "age_co": age_co,
        "tau_death_I1": tau_death_I1,
        "k_death_I1": k_death_I1,
        "tau_death_I2": tau_death_I2,
        "k_death_I2": k_death_I2,
        "tau_death_Co": tau_death_Co,
        "k_death_Co": k_death_Co,
        "f_death_Co1": f_death_Co1,
        "f_death_Co2": f_death_Co2,
        "k_bind_I1": k_bind_I1,
        "k_bind_I2": k_bind_I2,
        "k_bind_Co1": k_bind_Co1,
        "k_bind_Co2": k_bind_Co2,
        "mu_T": mu_T,
        "mu_I1": mu_I1,
        "mu_I2": mu_I2,
        "mu_Co": mu_Co,
        "k_deathT": k_deathT,
        "k_int_I1": k_int_I1,
        "k_int_I2": k_int_I2,
        "k_int_Co1": k_int_Co1,
        "k_int_Co2": k_int_Co2,
        "eta_I1": eta_I1,
        "eta_I2": eta_I2,
        "eta_Co1": eta_Co1,
        "eta_Co2": eta_Co2,
        "f_repl_I1": f_repl_I1,
        "f_repl_I2": f_repl_I2,
        "f_repl_Co1": f_repl_Co1,
        "f_repl_Co2": f_repl_Co2,
        "k_repl_I1": k_repl_I1,
        "k_repl_I2": k_repl_I2,
        "k_repl_Co1": k_repl_Co1,
        "k_repl_Co2": k_repl_Co2,
        "f_rel_I1": f_rel_I1,
        "f_rel_I2": f_rel_I2,
        "f_rel_Co": f_rel_Co,
        "k_rel_I1": k_rel_I1,
        "k_rel_I2": k_rel_I2,
        "k_rel_Co": k_rel_Co,
        "k_degrBV_1": k_d_V1,
        "k_degrBV_2": k_d_V2,
        "k_lys": k_lys,
        "randV2_I1": randV2_I1,
        "randV1_I2": randV1_I2,
        "k_d_N1": k_d_N1,
        "k_d_N2": k_d_N2,
        "K_s": K_s,
        "Ys_T": Ys_T,
        "Ys_I1": Ys_I1,
        "Ys_I2": Ys_I2,
        "Ys_C": Ys_C,
        "threshold": threshold,
        "rtol": rtol,
        "Dtau": Dtau,
        "hmax0": hmax0
    }
    
    return p