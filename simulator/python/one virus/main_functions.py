import numpy as np
import warnings

# @profile
def main_RK23(t_vect, x0, x0_bind, x0_nucl, D, r_bleed, Cin, Sin, p, sum_h, scaling):
    t0 = t_vect[0]
    t_final = t_vect[1]
    t = t0

    x = np.array(x0)
    x_bind = np.array(x0_bind)
    x_nucl = np.array(x0_nucl)

    # Extract parameters from parameter object 'p'
    startI1 = p['startI1']
    endI1 = p['endI1']
    startB1 = p['startB1']
    endB1 = p['endB1']
    n_bins = p['n_bins']
    n_bins_inf = p['n_bins_inf']
    age_bins = p['age_bins']
    tau_death_I1 = p['tau_death_I1']
    k_death_I1 = p['k_death_I1']
    k_bind_I1 = p['k_bind_I1']
    mu_T = p['mu_T']
    mu_I1 = p['mu_I1']
    k_deathT = p['k_deathT']
    k_int_I1 = p['k_int_I1']
    eta_I1 = p['eta_I1']
    f_repl_I1 = p['f_repl_I1']
    k_repl_I1 = p['k_repl_I1']
    f_rel_I1 = p['f_rel_I1']
    k_rel_I1 = p['k_rel_I1']
    k_degrBV_1 = p['k_degrBV_1']
    k_d_N1 = p['k_d_N1']
    k_lys = p['k_lys']
    K_s = p['K_s']
    Ys_T = p['Ys_T']
    Ys_I1 = p['Ys_I1']
    threshold = p['threshold']
    rtol = p['rtol']
    Dtau = p['Dtau']
    hmax0 = p['hmax0']


    # First step selection
    k1, k1_bind, k1_nucl = inf_model(t0, x0, x0_bind, x0_nucl, r_bleed, Cin, Sin, D, k_bind_I1, 
                  mu_T, mu_I1, k_deathT, k_death_I1, tau_death_I1, 
                  k_int_I1, eta_I1, k_repl_I1, f_repl_I1, f_rel_I1, 
                  k_rel_I1, k_degrBV_1, k_d_N1, K_s, Ys_T, Ys_I1, 
                  startI1, endI1, startB1, n_bins_inf, age_bins, 
                  n_bins, scaling, k_lys)

    r = max(np.linalg.norm(k1 / np.maximum(np.abs(x0), threshold), ord=np.inf),
            np.linalg.norm(k1_bind / np.maximum(np.abs(x0_bind), threshold), ord=np.inf))

    h = 0.8 * (rtol ** (1/3)) / r

    nofailed = True

    # Start loop
    while t < t_final:
        # Step control
        hmin = np.finfo(float).eps
        hmax = min(min(hmax0, Dtau - sum_h), t_final - t)
        h = max(hmin, min(h, hmax))

        # Intermediate steps (k2 and k3)
        k2, k2_bind, k2_nucl = inf_model(t + h/2, x + k1*h/2, x_bind + k1_bind*h/2, x_nucl + k1_nucl*h/2,
                                        r_bleed, Cin, Sin, D, k_bind_I1, 
                                        mu_T, mu_I1, k_deathT, k_death_I1, tau_death_I1, 
                                        k_int_I1, eta_I1, k_repl_I1, f_repl_I1, f_rel_I1, 
                                        k_rel_I1, k_degrBV_1, k_d_N1, K_s, Ys_T, Ys_I1, 
                                        startI1, endI1, startB1, n_bins_inf, age_bins, 
                                        n_bins, scaling, k_lys)

        k3, k3_bind, k3_nucl = inf_model(t + 3*h/4, x + 3*k2*h/4, x_bind + 3*k2_bind*h/4, x_nucl + 3*k2_nucl*h/4,
                                        r_bleed, Cin, Sin, D, k_bind_I1, 
                                        mu_T, mu_I1, k_deathT, k_death_I1, tau_death_I1, 
                                        k_int_I1, eta_I1, k_repl_I1, f_repl_I1, f_rel_I1, 
                                        k_rel_I1, k_degrBV_1, k_d_N1, K_s, Ys_T, Ys_I1, 
                                        startI1, endI1, startB1, n_bins_inf, age_bins, 
                                        n_bins, scaling, k_lys)

        # Update time and state variables
        t_new = t + h
        x_new = np.real(x + h * (2*k1 + 3*k2 + 4*k3) / 9)
        x_bind_new = np.real(x_bind + h * (2*k1_bind + 3*k2_bind + 4*k3_bind) / 9)
        x_nucl_new = np.real(x_nucl + h * (2*k1_nucl + 3*k2_nucl + 4*k3_nucl) / 9)

        # Final step (k4)
        k4, k4_bind, k4_nucl = inf_model(t_new, x_new, x_bind_new, x_nucl_new,
                                         r_bleed, Cin, Sin, D, k_bind_I1, 
                                         mu_T, mu_I1, k_deathT, k_death_I1, tau_death_I1, 
                                         k_int_I1, eta_I1, k_repl_I1, f_repl_I1, f_rel_I1, 
                                         k_rel_I1, k_degrBV_1, k_d_N1, K_s, Ys_T, Ys_I1, 
                                         startI1, endI1, startB1, n_bins_inf, age_bins, 
                                         n_bins, scaling, k_lys)

        # Estimate error
        e = h * (-5*k1 + 6*k2 + 8*k3 - 9*k4) / 72
        e_bind = h * (-5*k1_bind + 6*k2_bind + 8*k3_bind - 9*k4_bind) / 72
        e_nucl = h * (-5*k1_nucl + 6*k2_nucl + 8*k3_nucl - 9*k4_nucl) / 72

        err_est = max(np.linalg.norm(e / np.maximum(np.maximum(np.abs(x), np.abs(x_new)), threshold), ord=np.inf),
                      np.linalg.norm(e_bind / np.maximum(np.maximum(np.abs(x_bind), np.abs(x_bind_new)), threshold), ord=np.inf),
                      np.linalg.norm(e_nucl / np.maximum(np.maximum(np.abs(x_nucl), np.abs(x_nucl_new)), threshold), ord=np.inf))

        if err_est <= rtol:
            # Accept step
            t = t_new
            x = x_new
            x_bind = x_bind_new
            x_nucl = x_nucl_new

            x[x < 0] = 0
            x_bind[x_bind < 0] = 0
            x_nucl[x_nucl < 0] = 0
        
            sum_h += h
        
            k1 = k4
            k1_bind = k4_bind
            k1_nucl = k4_nucl
        
            nofailed = True
            
            if sum_h >= Dtau or t >= t_final:
                t = round(t, 2)
                
                sum_h = 0
            
                # Reallocate I1
                x[endI1] += x[endI1-1]
                x[startI1+1:endI1] = x[startI1:endI1-1]
                x[startI1] = 0
                
                # Reallocate virus bound to I1
                x_bind[endB1] += x_bind[endB1-1]
                x_bind[startB1+1:endB1] = x_bind[startB1:endB1-1]
                x_bind[startB1] = 0
            
                # Reallocate nuclear virus in I1
                x_nucl[endB1] += x_nucl[endB1-1]
                x_nucl[startB1+1:endB1] = x_nucl[startB1:endB1-1]
                x_nucl[startB1] = 0
            
                # Update first slope after re-allocation
                k1, k1_bind, k1_nucl = inf_model(
                    t, x, x_bind, x_nucl,
                    r_bleed, Cin, Sin, D, k_bind_I1, 
                    mu_T, mu_I1, k_deathT, k_death_I1, tau_death_I1, 
                    k_int_I1, eta_I1, k_repl_I1, f_repl_I1, f_rel_I1, 
                    k_rel_I1, k_degrBV_1, k_d_N1, K_s, Ys_T, Ys_I1, 
                    startI1, endI1, startB1, n_bins_inf, age_bins, 
                    n_bins, scaling, k_lys
                )

        else:
            # Reject step and reduce step size
            nofailed = False
        
        # Adjust the next step size
        if nofailed:
            nofailed = False
            h = max(hmin, h * max(0.5, 0.8 * (rtol / err_est)**(1/3)))
        else:
            h = max(hmin, 0.5 * h)
        # Avoid step size getting too small
        if h <= hmin:
            warnings.warn('Step size became too small. Integration may not be accurate.', RuntimeWarning)
            break
        
    return t,x,x_bind,x_nucl,sum_h

def inf_model(t, x, x_bind, x_nucl, r_bleed, Cin, Sin, D, k_bind_I1, 
                        mu_T, mu_I1, k_deathT, k_death_I1, tau_death_I1, 
                        k_int_I1, eta_I1, k_repl_I1, f_repl_I1, f_rel_I1, 
                        k_rel_I1, k_degrBV_1, k_d_N1, K_s, Ys_T, Ys_I1, 
                        startI1, endI1, startB1, n_bins_inf, age_bins, 
                        n_bins, scaling, k_lys):
    
    dxdt = np.zeros_like(x)  
    b = r_bleed * D
    
    # Species extraction
    T = x[0]
    V1 = x[1]
    I1 = x[startI1:endI1 + 1]
    NV = x[-2]
    S = x[-1]

    # Pre-allocate dx_bind_dt and dx_nucl_dt
    dx_bind_dt = np.zeros_like(x_bind)
    dx_nucl_dt = np.zeros_like(x_nucl)

    totI1 = np.sum(I1)  
    rel_v1 = 0  # Release scalar

    growth_lim = max(0, S / (K_s + S))
    growthT = mu_T * T * growth_lim
    dxdt[0] = growthT - k_deathT * T - k_bind_I1[0] * T * V1 + Cin * D / scaling - b * T
    dxdt[-2] += k_deathT * T
    dxdt[startI1] += k_bind_I1[0] * T * V1
    dx_bind_dt[startB1] += k_bind_I1[0] * T * V1

    ages = age_bins  
    growthI1 = mu_I1 * I1 * growth_lim
    dxdt[startI1:startI1 + n_bins] += growthI1

    # Calculate death rates based on DNA copy numbers and age
    death_age_factor = (0.5 + 0.5 * np.tanh((ages - tau_death_I1) / 0.3))
    k_deathI1_vector = k_deathT * (0.5 + 0.5 * np.tanh((tau_death_I1 - ages) / 0.3)) + \
                       k_death_I1 * np.maximum(1, np.log((x_nucl[:n_bins]) / (I1 + 1e-50))) * death_age_factor

    deathI1 = k_deathI1_vector * I1
    dxdt[startI1:startI1 + n_bins] -= deathI1
    dxdt[-2] += np.sum(deathI1)  # Sum of death rates
    
    # Virus uptake
    I1_uptake = k_bind_I1[:n_bins_inf] * I1[:n_bins_inf] * V1
    dx_bind_dt[:n_bins_inf] += I1_uptake
    bind1 = np.sum(I1_uptake) + k_bind_I1[0] * T * V1

    # Virus internalization and replication 
    Ein1 = k_int_I1 * x_bind
    dx_bind_dt[:n_bins] -= Ein1 + k_deathI1_vector * x_bind
    dx_nucl_dt[:n_bins] += Ein1 * eta_I1 - k_deathI1_vector * x_nucl[:n_bins] - k_d_N1 * x_nucl[:n_bins]
    
    # DNA replication and viral release
    dx_nucl_dt[:n_bins] += k_repl_I1 * x_nucl[:n_bins] * f_repl_I1
    rel_v1 = np.sum(k_rel_I1 * I1 * f_rel_I1)
    
    # Virions balance
    dxdt[1] = -bind1 - k_degrBV_1 * V1 + rel_v1 - D * V1

    # Cells: bleeding and outlet 
    dxdt[startI1:startI1 + n_bins] -= b * I1
    dx_bind_dt -= b * x_bind
    dx_nucl_dt -= b * x_nucl

    # Nonviable cells and substrate consumption
    dxdt[-2] -= (k_lys + b) * NV
    dx_nucl_dt[-1] = np.sum(k_deathI1_vector*x_nucl[:n_bins]) - (k_lys + b + k_d_N1) * x_nucl[-1]
    dxdt[-1] = D * (Sin - S) - (Ys_T * T + Ys_I1 * totI1) * S / (1e-2 / scaling + S)

    return dxdt, dx_bind_dt, dx_nucl_dt

