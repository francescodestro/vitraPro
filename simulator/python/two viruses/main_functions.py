import numpy as np

# @profile
def main_RK23(t_vect, x0, x0_bind, x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling):
    
    t0 = t_vect[0]
    t_final = t_vect[1]
    t = t0
    
    x = x0.copy()
    x_bind = x0_bind.copy()
    x_nucl = x0_nucl.copy()
    
    # Get parameters from parameters dictionary
    startI1 = p['startI1']
    endI1 = p['endI1']
    startI2 = p['startI2']
    endI2 = p['endI2']
    startCo = p['startCo']
    endCo = p['endCo']
    startB1 = p['startB1']
    startB2 = p['startB2']
    endB1 = p['endB1']
    endB2 = p['endB2']
    startB1Co = p['startB1Co']
    startB2Co = p['startB2Co']
    ind_scCo = p['ind_scCo']
    age_end_inf = p['age_end_inf']
    a = p['a']
    
    ind = p['ind']
    ind_scB1Co = p['ind_scB1Co']
    ind_scB2Co = p['ind_scB2Co']
    
    n_bins = p['n_bins']
    
    n_bins_inf = p['n_bins_inf']
    age_bins = p['age_bins']
    age_co = p['age_co']
    
    tau_death_I1 = p['tau_death_I1']
    k_death_I1 = p['k_death_I1']
    tau_death_I2 = p['tau_death_I2']
    k_death_I2 = p['k_death_I2']
    tau_death_Co = p['tau_death_Co']
    k_death_Co = p['k_death_Co']
    f_death_Co1 = p['f_death_Co1']
    f_death_Co2 = p['f_death_Co2']
    
    k_bind_I1 = p['k_bind_I1']
    k_bind_I2 = p['k_bind_I2']
    k_bind_Co1 = p['k_bind_Co1']
    k_bind_Co2 = p['k_bind_Co2']
    
    mu_T = p['mu_T']
    mu_I1 = p['mu_I1']
    mu_I2 = p['mu_I2']
    mu_Co = p['mu_Co']
    
    k_deathT = p['k_deathT']
    
    k_int_I1 = p['k_int_I1']
    k_int_I2 = p['k_int_I2']
    k_int_Co1 = p['k_int_Co1']
    k_int_Co2 = p['k_int_Co2']
    
    eta_I1 = p['eta_I1']
    eta_I2 = p['eta_I2']
    eta_Co1 = p['eta_Co1']
    eta_Co2 = p['eta_Co2']
    
    f_repl_I1 = p['f_repl_I1']
    f_repl_I2 = p['f_repl_I2']
    f_repl_Co1 = p['f_repl_Co1']
    f_repl_Co2 = p['f_repl_Co2']
    
    k_repl_I1 = p['k_repl_I1']
    k_repl_I2 = p['k_repl_I2']
    k_repl_Co1 = p['k_repl_Co1']
    k_repl_Co2 = p['k_repl_Co2']
    
    f_rel_I1 = p['f_rel_I1']
    f_rel_I2 = p['f_rel_I2']
    f_rel_Co = p['f_rel_Co']
    k_rel_I1 = p['k_rel_I1']
    k_rel_I2 = p['k_rel_I2']
    k_rel_Co = p['k_rel_Co']
    
    k_degrBV_1 = p['k_degrBV_1']
    k_degrBV_2 = p['k_degrBV_2']
    
    randV2_I1 = p['randV2_I1']
    randV1_I2 = p['randV1_I2']
    
    k_d_N1 = p['k_d_N1']
    k_d_N2 = p['k_d_N2']
    
    k_lys = p['k_lys']
    
    K_s = p['K_s']
    Ys_T = p['Ys_T']
    Ys_I1 = p['Ys_I1']
    Ys_I2 = p['Ys_I2']
    Ys_C = p['Ys_C']
    
    threshold = p['threshold']
    rtol = p['rtol']
    Dtau = p['Dtau']
    hmax0 = p['hmax0']
    
    # first step selection
    k1, k1_bind, k1_nucl = inf_model(t0, x0, x0_bind, x0_nucl, r_bleed, Cin, 
         Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, mu_I1, 
         mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, f_death_Co1, 
         f_death_Co2, tau_death_I1, tau_death_I2, tau_death_Co, k_int_I1, 
         k_int_I2, k_int_Co1, eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, 
         k_repl_I1, k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, f_repl_I2, 
         f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, f_rel_Co, k_rel_I1, 
         k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1, k_d_N2, K_s, 
         Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, randV1_I2, startI1, endI1, 
         startI2, endI2, startCo, endCo, startB1, startB2, n_bins_inf, 
         age_bins, ind_scB1Co, ind_scB2Co, n_bins, age_co, ind, startB1Co, 
         startB2Co, scaling, k_lys)
    r = max(np.linalg.norm(k1 / np.maximum(np.abs(x0), threshold), np.inf), 
        np.linalg.norm(k1_bind / np.maximum(np.abs(x0_bind), threshold), np.inf))   
    h = 0.8 * (rtol ** (1/3)) / r

    nofailed = True
    
    while t < round(t_final, 2):
        # Step control
        hmin = np.finfo(float).eps
        hmax = min(min(hmax0, Dtau - sum_h), t_final - t)
        h = max(hmin, min(h, hmax))
        
        # Calculate k2, k3
        k2, k2_bind, k2_nucl = inf_model(t + h/2, x + k1*h/2, x_bind + k1_bind*h/2, 
             x_nucl + k1_nucl*h/2, r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, 
             k_bind_Co1, k_bind_Co2, mu_T, mu_I1, mu_I2, mu_Co, k_deathT, 
             k_death_I1, k_death_I2, k_death_Co, f_death_Co1, f_death_Co2, 
             tau_death_I1, tau_death_I2, tau_death_Co, k_int_I1, k_int_I2, 
             k_int_Co1, eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, 
             k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, 
             f_repl_Co2, f_rel_I1, f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, 
             k_degrBV_1, k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, 
             randV2_I1, randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, 
             startB1, startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)
        
        k3, k3_bind, k3_nucl = inf_model(t + 3*h/4, x + 3*k2*h/4, x_bind + 3*k2_bind*h/4,
             x_nucl + 3*k2_nucl*h/4, r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, 
             k_bind_Co1, k_bind_Co2, mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1,
             k_death_I2, k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, 
             tau_death_I2, tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, 
             eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, 
             k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, randV1_I2, 
             startI1, endI1, startI2, endI2, startCo, endCo, startB1, startB2, 
             n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, age_co, ind, 
             startB1Co, startB2Co, scaling, k_lys)
        
        t_new = t + h
 
        x_new = np.real(x + h * (2*k1 + 3*k2 + 4*k3) / 9)
        x_bind_new = np.real(x_bind + h * (2*k1_bind + 3*k2_bind + 4*k3_bind) / 9)
        x_nucl_new = np.real(x_nucl + h * (2*k1_nucl + 3*k2_nucl + 4*k3_nucl) / 9)
        
        k4, k4_bind, k4_nucl = inf_model(t_new, x_new, x_bind_new, x_nucl_new, 
             r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
             mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
             k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, 
             tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, eta_I2, 
             eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, 
             k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
             randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, startB1, 
             startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)

        # Estimate error
        e = h * (-5*k1 + 6*k2 + 8*k3 - 9*k4) / 72
        e_bind = h * (-5*k1_bind + 6*k2_bind + 8*k3_bind - 9*k4_bind) / 72
        e_nucl = h * (-5*k1_nucl + 6*k2_nucl + 8*k3_nucl - 9*k4_nucl) / 72
        
        err_est = max(
            np.linalg.norm(e / np.maximum(np.maximum(np.abs(x), np.abs(x_new)), threshold), np.inf),
            np.linalg.norm(e_bind / np.maximum(np.maximum(np.abs(x_bind), np.abs(x_bind_new)), threshold), np.inf),
            np.linalg.norm(e_nucl / np.maximum(np.maximum(np.abs(x_nucl), np.abs(x_nucl_new)), threshold), np.inf)
        )
        
        if err_est <= rtol:
            t = t_new
            x = x_new.copy()
            x_bind = x_bind_new.copy()
            x_nucl = x_nucl_new.copy()
        
            x[x < 0] = 0
            x_bind[x_bind < 0] = 0
            x_nucl[x_nucl < 0] = 0
        
            k1 = k4
            k1_bind = k4_bind
            k1_nucl = k4_nucl
        
            sum_h = sum_h + h
        
            temp = 1.25 * ((err_est / rtol) ** (1/3))
            h = h * min(5, 1 / temp)
            nofailed = True
            
            # Age switch
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
            
                # Reallocate I2
                x[endI2] += x[endI2-1]
                x[startI2+1:endI2] = x[startI2:endI2-1]
                x[startI2] = 0
            
                # Reallocate virus bound to I2
                x_bind[endB2] += x_bind[endB2-1]
                x_bind[startB2+1:endB2] = x_bind[startB2:endB2-1]
                x_bind[startB2] = 0
            
                # Reallocate nuclear virus in I2
                x_nucl[endB2] += x_nucl[endB2-1]
                x_nucl[startB2+1:endB2] = x_nucl[startB2:endB2-1]
                x_nucl[startB2] = 0
            
                # Reallocate co-infected cells and virus bound to them
                # Max age bins
                i1 = n_bins-1
                i2=np.arange(max(1, n_bins - n_bins_inf), n_bins)
                x[ind_scCo[i2+n_bins*i1]] += x[ind_scCo[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB1Co[i2+n_bins*i1]] += x_bind[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB2Co[i2+n_bins*i1]] += x_bind[ind_scB2Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB1Co[i2+n_bins*i1]] += x_nucl[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB2Co[i2+n_bins*i1]] += x_nucl[ind_scB2Co[i2-1+n_bins*(i1-1)]]
            
                i1 = np.arange(max(1, n_bins - n_bins_inf), n_bins - 1)
                i2 = n_bins-1
                x[ind_scCo[i2+n_bins*i1]] += x[ind_scCo[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB1Co[i2+n_bins*i1]] += x_bind[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB2Co[i2+n_bins*i1]] += x_bind[ind_scB2Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB1Co[i2+n_bins*i1]] += x_nucl[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB2Co[i2+n_bins*i1]] += x_nucl[ind_scB2Co[i2-1+n_bins*(i1-1)]]
            
                # Translation of intermediate bins
                for i2 in np.arange(n_bins - 2, n_bins_inf-1, -1):
                    for i1 in np.arange(i2, i2 - n_bins_inf, -1):
                        x[ind_scCo[i2 + n_bins * i1]] = x[ind_scCo[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB1Co[i2 + n_bins * i1]] = x_bind[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB2Co[i2 + n_bins * i1]] = x_bind[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB1Co[i2 + n_bins * i1]] = x_nucl[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB2Co[i2 + n_bins * i1]] = x_nucl[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                
                for i1 in np.arange(n_bins - 2, n_bins_inf-1, -1):
                    for i2 in np.arange(i1 - 1, i1 - n_bins_inf, -1):
                        x[ind_scCo[i2 + n_bins * i1]] = x[ind_scCo[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB1Co[i2 + n_bins * i1]] = x_bind[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB2Co[i2 + n_bins * i1]] = x_bind[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB1Co[i2 + n_bins * i1]] = x_nucl[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB2Co[i2 + n_bins * i1]] = x_nucl[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                
                for i1 in np.arange(n_bins_inf + a - 1, 0, -1):
                    for i2 in np.arange(n_bins_inf + a - 1, 0, -1):
                        x[ind_scCo[i2 + n_bins * i1]] = x[ind_scCo[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB1Co[i2 + n_bins * i1]] = x_bind[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB2Co[i2 + n_bins * i1]] = x_bind[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB1Co[i2 + n_bins * i1]] = x_nucl[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB2Co[i2 + n_bins * i1]] = x_nucl[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        
                        
                # Reset first bin
                i1 = 0
                i2 = np.arange(0, n_bins_inf) 
                x[ind_scCo[i2 + n_bins * i1]] = 0
                x_bind[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_bind[ind_scB2Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB2Co[i2 + n_bins * i1]] = 0
                
                i1 = np.arange(1, n_bins_inf)
                i2 = 0
                x[ind_scCo[i2 + n_bins * i1]] = 0
                x_bind[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_bind[ind_scB2Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB2Co[i2 + n_bins * i1]] = 0
                
                # Update first slope after re-allocation
                k1, k1_bind, k1_nucl = inf_model(t, x, x_bind, x_nucl, r_bleed, 
                    Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
                    mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
                    k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, 
                    tau_death_I2, tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, 
                    eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, 
                    k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, f_repl_I2, 
                    f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, f_rel_Co, 
                    k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, 
                    k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
                    randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, 
                    startB1, startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, 
                    n_bins, age_co, ind, startB1Co, startB2Co, scaling, k_lys)
                
        else:
            if nofailed:
                nofailed = False
                h = max(hmin, h * max(0.5, 0.8*(rtol/err_est)**(1/3)))
            else:
                h = max(hmin, 0.5 * h)
        
        if h <= hmin:
            print('Step size {} too small at t = {}\n'.format(h, t))
            break
     
    return  t,x,x_bind,x_nucl,sum_h

def main_RK45(t_vect, x0, x0_bind, x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling):
    
    t0 = t_vect[0]
    t_final = t_vect[1]
    t = t0
    
    x = x0.copy()
    x_bind = x0_bind.copy()
    x_nucl = x0_nucl.copy()
    
    # RK45 parameters
    a2 = 1/5
    a3 = 3/10
    a4 = 4/5
    a5 = 8/9
    
    b11 = 1/5
    b21 = 3/40
    b31 = 44/45
    b41 = 19372/6561
    b51 = 9017/3168
    b61 = 35/384
    b22 = 9/40
    b32 = -56/15
    b42 = -25360/2187
    b52 = -355/33
    b33 = 32/9
    b43 = 64448/6561
    b53 = 46732/5247
    b63 = 500/1113
    b44 = -212/729
    b54 = 49/176
    b64 = 125/192
    b55 = -5103/18656
    b65 = -2187/6784
    b66 = 11/84
    
    e1 = 71/57600
    e3 = -71/16695
    e4 = 71/1920
    e5 = -17253/339200
    e6 = 22/525
    e7 = -1/40
    
    # Get parameters from parameters dictionary
    startI1 = p['startI1']
    endI1 = p['endI1']
    startI2 = p['startI2']
    endI2 = p['endI2']
    startCo = p['startCo']
    endCo = p['endCo']
    startB1 = p['startB1']
    startB2 = p['startB2']
    endB1 = p['endB1']
    endB2 = p['endB2']
    startB1Co = p['startB1Co']
    startB2Co = p['startB2Co']
    ind_scCo = p['ind_scCo']
    age_end_inf = p['age_end_inf']
    a = p['a']
    
    ind = p['ind']
    ind_scB1Co = p['ind_scB1Co']
    ind_scB2Co = p['ind_scB2Co']
    
    n_bins = p['n_bins']
    
    n_bins_inf = p['n_bins_inf']
    age_bins = p['age_bins']
    age_co = p['age_co']
    
    tau_death_I1 = p['tau_death_I1']
    k_death_I1 = p['k_death_I1']
    tau_death_I2 = p['tau_death_I2']
    k_death_I2 = p['k_death_I2']
    tau_death_Co = p['tau_death_Co']
    k_death_Co = p['k_death_Co']
    f_death_Co1 = p['f_death_Co1']
    f_death_Co2 = p['f_death_Co2']
    
    k_bind_I1 = p['k_bind_I1']
    k_bind_I2 = p['k_bind_I2']
    k_bind_Co1 = p['k_bind_Co1']
    k_bind_Co2 = p['k_bind_Co2']
    
    mu_T = p['mu_T']
    mu_I1 = p['mu_I1']
    mu_I2 = p['mu_I2']
    mu_Co = p['mu_Co']
    
    k_deathT = p['k_deathT']
    
    k_int_I1 = p['k_int_I1']
    k_int_I2 = p['k_int_I2']
    k_int_Co1 = p['k_int_Co1']
    k_int_Co2 = p['k_int_Co2']
    
    eta_I1 = p['eta_I1']
    eta_I2 = p['eta_I2']
    eta_Co1 = p['eta_Co1']
    eta_Co2 = p['eta_Co2']
    
    f_repl_I1 = p['f_repl_I1']
    f_repl_I2 = p['f_repl_I2']
    f_repl_Co1 = p['f_repl_Co1']
    f_repl_Co2 = p['f_repl_Co2']
    
    k_repl_I1 = p['k_repl_I1']
    k_repl_I2 = p['k_repl_I2']
    k_repl_Co1 = p['k_repl_Co1']
    k_repl_Co2 = p['k_repl_Co2']
    
    f_rel_I1 = p['f_rel_I1']
    f_rel_I2 = p['f_rel_I2']
    f_rel_Co = p['f_rel_Co']
    k_rel_I1 = p['k_rel_I1']
    k_rel_I2 = p['k_rel_I2']
    k_rel_Co = p['k_rel_Co']
    
    k_degrBV_1 = p['k_degrBV_1']
    k_degrBV_2 = p['k_degrBV_2']
    
    randV2_I1 = p['randV2_I1']
    randV1_I2 = p['randV1_I2']
    
    k_d_N1 = p['k_d_N1']
    k_d_N2 = p['k_d_N2']
    
    k_lys = p['k_lys']
    
    K_s = p['K_s']
    Ys_T = p['Ys_T']
    Ys_I1 = p['Ys_I1']
    Ys_I2 = p['Ys_I2']
    Ys_C = p['Ys_C']
    
    threshold = p['threshold']
    rtol = p['rtol']
    Dtau = p['Dtau']
    hmax0 = p['hmax0']
    
    # first step selection
    k1, k1_bind, k1_nucl = inf_model(t0, x0, x0_bind, x0_nucl, r_bleed, Cin, 
         Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, mu_I1, 
         mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, f_death_Co1, 
         f_death_Co2, tau_death_I1, tau_death_I2, tau_death_Co, k_int_I1, 
         k_int_I2, k_int_Co1, eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, 
         k_repl_I1, k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, f_repl_I2, 
         f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, f_rel_Co, k_rel_I1, 
         k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1, k_d_N2, K_s, 
         Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, randV1_I2, startI1, endI1, 
         startI2, endI2, startCo, endCo, startB1, startB2, n_bins_inf, 
         age_bins, ind_scB1Co, ind_scB2Co, n_bins, age_co, ind, startB1Co, 
         startB2Co, scaling, k_lys)
    r = max(np.linalg.norm(k1 / np.maximum(np.abs(x0), threshold), np.inf), 
        np.linalg.norm(k1_bind / np.maximum(np.abs(x0_bind), threshold), np.inf))   
    h = 0.8 * (rtol ** (1/3)) / r

    nofailed = True
    
    while t < round(t_final, 2):
        # Step control
        hmin = np.finfo(float).eps
        hmax = min(min(hmax0, Dtau - sum_h), t_final - t)
        h = max(hmin, min(h, hmax))
        
        # step 2
        x2 = x + h * (b11*k1 )
        x_bind2 = x_bind + h * (b11*k1_bind )
        x_nucl2 = x_nucl + h * (b11*k1_nucl )
        t2 = t + h * a2
        
        # Calculate k2, k3
        k2, k2_bind, k2_nucl = inf_model(t2, x2, x_bind2, x_nucl2, 
             r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, 
             k_bind_Co1, k_bind_Co2, mu_T, mu_I1, mu_I2, mu_Co, k_deathT, 
             k_death_I1, k_death_I2, k_death_Co, f_death_Co1, f_death_Co2, 
             tau_death_I1, tau_death_I2, tau_death_Co, k_int_I1, k_int_I2, 
             k_int_Co1, eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, 
             k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, 
             f_repl_Co2, f_rel_I1, f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, 
             k_degrBV_1, k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, 
             randV2_I1, randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, 
             startB1, startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)
        
        # step3
        x3 = x + h * (b21*k1 + b22*k2 )
        x_bind3 = x_bind + h * (b21*k1_bind + b22*k2_bind)
        x_nucl3 = x_nucl + h * (b21*k1_nucl + b22*k2_nucl)
        t3 = t + h * a3
        
        k3, k3_bind, k3_nucl = inf_model(t3, x3, x_bind3,
             x_nucl3, r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, 
             k_bind_Co1, k_bind_Co2, mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1,
             k_death_I2, k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, 
             tau_death_I2, tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, 
             eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, 
             k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, randV1_I2, 
             startI1, endI1, startI2, endI2, startCo, endCo, startB1, startB2, 
             n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, age_co, ind, 
             startB1Co, startB2Co, scaling, k_lys)
        
        # step4
        x4 = x + h * (b31*k1 + b32*k2 + b33*k3)
        x_bind4 = x_bind + h * (b31*k1_bind + b32*k2_bind + b33*k3_bind)
        x_nucl4 = x_nucl + h * (b31*k1_nucl + b32*k2_nucl + b33*k3_nucl)
        t4 = t + h * a4
        
        k4, k4_bind, k4_nucl = inf_model(t4, x4, x_bind4, x_nucl4, 
             r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
             mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
             k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, 
             tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, eta_I2, 
             eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, 
             k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
             randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, startB1, 
             startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)
        
        # step5
        x5 = x + h * (b41*k1 + b42*k2 + b43*k3 + b44*k4 )
        x_bind5 = x_bind + h * (b41*k1_bind + b42*k2_bind + 
            b43*k3_bind + b44*k4_bind )
        x_nucl5 = x_nucl + h * (b41*k1_nucl + b42*k2_nucl + 
            b43*k3_nucl + b44*k4_nucl )
        t5 = t + h * a5
        
        k5, k5_bind, k5_nucl = inf_model(t5, x5, x_bind5, x_nucl5, 
             r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
             mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
             k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, 
             tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, eta_I2, 
             eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, 
             k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
             randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, startB1, 
             startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)
        
        # step6
        x6 = x + h * (b51*k1 + b52*k2 + b53*k3 + b54*k4 + b55*k5 )
        x_bind6 = x_bind + h * (b51*k1_bind + b52*k2_bind + b53*k3_bind + 
            b54*k4_bind + b55*k5_bind )
        x_nucl6 = x_nucl + h * (b51*k1_nucl + b52*k2_nucl + b53*k3_nucl + 
            b54*k4_nucl + b55*k5_nucl )
        t6 = t + h
        
        k6, k6_bind, k6_nucl = inf_model(t6, x6, x_bind6, x_nucl6, 
             r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
             mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
             k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, 
             tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, eta_I2, 
             eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, 
             k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
             randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, startB1, 
             startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)
        
        # Final step
        t_new=t+h
        
        x_new=np.real(x + h* (b61*k1 + b63*k3 + b64*k4 + b65*k5 + b66*k6))
        x_bind_new=np.real(x_bind + h* (b61*k1_bind + b63*k3_bind + 
            b64*k4_bind + b65*k5_bind + b66*k6_bind))
        x_nucl_new=np.real(x_nucl + h* (b61*k1_nucl + b63*k3_nucl + 
            b64*k4_nucl + b65*k5_nucl + b66*k6_nucl))
        
        k7, k7_bind, k7_nucl = inf_model(t_new, x_new, x_bind_new, x_nucl_new, 
             r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
             mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
             k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, 
             tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, eta_I1, eta_I2, 
             eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, k_repl_I2, k_repl_Co1, 
             k_repl_Co2, f_repl_I1, f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, 
             f_rel_I2, f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, 
             k_degrBV_2, k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
             randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, startB1, 
             startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, n_bins, 
             age_co, ind, startB1Co, startB2Co, scaling, k_lys)

        # estimate error
        e = h* (k1*e1 + k3*e3 + k4*e4 + k5*e5 + k6*e6 + k7*e7)
        e_bind = h* (k1_bind*e1 + k3_bind*e3 + k4_bind*e4 + k5_bind*e5 + 
            k6_bind*e6 + k7_bind*e7)
        e_nucl = h* (k1_nucl*e1 + k3_nucl*e3 + k4_nucl*e4 + k5_nucl*e5 + 
            k6_nucl*e6 + k7_nucl*e7)
    
        
        err_est = max(
            np.linalg.norm(e / np.maximum(np.maximum(np.abs(x), np.abs(x_new)), threshold), np.inf),
            np.linalg.norm(e_bind / np.maximum(np.maximum(np.abs(x_bind), np.abs(x_bind_new)), threshold), np.inf),
            np.linalg.norm(e_nucl / np.maximum(np.maximum(np.abs(x_nucl), np.abs(x_nucl_new)), threshold), np.inf)
        )
        
        if err_est <= rtol:
            t = t_new
            x = x_new.copy()
            x_bind = x_bind_new.copy()
            x_nucl = x_nucl_new.copy()
        
            x[x < 0] = 0
            x_bind[x_bind < 0] = 0
            x_nucl[x_nucl < 0] = 0
        
            k1 = k4
            k1_bind = k4_bind
            k1_nucl = k4_nucl
        
            sum_h = sum_h + h
        
            temp = 1.25 * ((err_est / rtol) ** (1/3))
            h = h * min(5, 1 / temp)
            nofailed = True
            
            # Age switch
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
            
                # Reallocate I2
                x[endI2] += x[endI2-1]
                x[startI2+1:endI2] = x[startI2:endI2-1]
                x[startI2] = 0
            
                # Reallocate virus bound to I2
                x_bind[endB2] += x_bind[endB2-1]
                x_bind[startB2+1:endB2] = x_bind[startB2:endB2-1]
                x_bind[startB2] = 0
            
                # Reallocate nuclear virus in I2
                x_nucl[endB2] += x_nucl[endB2-1]
                x_nucl[startB2+1:endB2] = x_nucl[startB2:endB2-1]
                x_nucl[startB2] = 0
            
                # Reallocate co-infected cells and virus bound to them
                # Max age bins
                i1 = n_bins-1
                i2=np.arange(max(1, n_bins - n_bins_inf), n_bins)
                x[ind_scCo[i2+n_bins*i1]] += x[ind_scCo[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB1Co[i2+n_bins*i1]] += x_bind[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB2Co[i2+n_bins*i1]] += x_bind[ind_scB2Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB1Co[i2+n_bins*i1]] += x_nucl[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB2Co[i2+n_bins*i1]] += x_nucl[ind_scB2Co[i2-1+n_bins*(i1-1)]]
            
                i1 = np.arange(max(1, n_bins - n_bins_inf), n_bins - 1)
                i2 = n_bins-1
                x[ind_scCo[i2+n_bins*i1]] += x[ind_scCo[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB1Co[i2+n_bins*i1]] += x_bind[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_bind[ind_scB2Co[i2+n_bins*i1]] += x_bind[ind_scB2Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB1Co[i2+n_bins*i1]] += x_nucl[ind_scB1Co[i2-1+n_bins*(i1-1)]]
                x_nucl[ind_scB2Co[i2+n_bins*i1]] += x_nucl[ind_scB2Co[i2-1+n_bins*(i1-1)]]
            
                # Translation of intermediate bins
                for i2 in np.arange(n_bins - 2, n_bins_inf-1, -1):
                    for i1 in np.arange(i2, i2 - n_bins_inf, -1):
                        x[ind_scCo[i2 + n_bins * i1]] = x[ind_scCo[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB1Co[i2 + n_bins * i1]] = x_bind[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB2Co[i2 + n_bins * i1]] = x_bind[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB1Co[i2 + n_bins * i1]] = x_nucl[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB2Co[i2 + n_bins * i1]] = x_nucl[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                
                for i1 in np.arange(n_bins - 2, n_bins_inf-1, -1):
                    for i2 in np.arange(i1 - 1, i1 - n_bins_inf, -1):
                        x[ind_scCo[i2 + n_bins * i1]] = x[ind_scCo[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB1Co[i2 + n_bins * i1]] = x_bind[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB2Co[i2 + n_bins * i1]] = x_bind[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB1Co[i2 + n_bins * i1]] = x_nucl[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB2Co[i2 + n_bins * i1]] = x_nucl[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                
                for i1 in np.arange(n_bins_inf + a - 1, 0, -1):
                    for i2 in np.arange(n_bins_inf + a - 1, 0, -1):
                        x[ind_scCo[i2 + n_bins * i1]] = x[ind_scCo[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB1Co[i2 + n_bins * i1]] = x_bind[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_bind[ind_scB2Co[i2 + n_bins * i1]] = x_bind[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB1Co[i2 + n_bins * i1]] = x_nucl[ind_scB1Co[i2 - 1 + n_bins * (i1 - 1)]]
                        x_nucl[ind_scB2Co[i2 + n_bins * i1]] = x_nucl[ind_scB2Co[i2 - 1 + n_bins * (i1 - 1)]]
                        
                        
                # Reset first bin
                i1 = 0
                i2 = np.arange(0, n_bins_inf) 
                x[ind_scCo[i2 + n_bins * i1]] = 0
                x_bind[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_bind[ind_scB2Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB2Co[i2 + n_bins * i1]] = 0
                
                i1 = np.arange(1, n_bins_inf)
                i2 = 0
                x[ind_scCo[i2 + n_bins * i1]] = 0
                x_bind[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_bind[ind_scB2Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB1Co[i2 + n_bins * i1]] = 0
                x_nucl[ind_scB2Co[i2 + n_bins * i1]] = 0
                
                # Update first slope after re-allocation
                k1, k1_bind, k1_nucl = inf_model(t, x, x_bind, x_nucl, r_bleed, 
                    Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, 
                    mu_T, mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, 
                    k_death_Co, f_death_Co1, f_death_Co2, tau_death_I1, 
                    tau_death_I2, tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, 
                    eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, 
                    k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, f_repl_I2, 
                    f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, f_rel_Co, 
                    k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, 
                    k_d_N1, k_d_N2, K_s, Ys_T, Ys_I1, Ys_I2, Ys_C, randV2_I1, 
                    randV1_I2, startI1, endI1, startI2, endI2, startCo, endCo, 
                    startB1, startB2, n_bins_inf, age_bins, ind_scB1Co, ind_scB2Co, 
                    n_bins, age_co, ind, startB1Co, startB2Co, scaling, k_lys)
                
        else:
            if nofailed:
                nofailed = False
                h = max(hmin, h * max(0.5, 0.8*(rtol/err_est)**(1/3)))
            else:
                h = max(hmin, 0.5 * h)
        
        if h <= hmin:
            print('Step size {} too small at t = {}\n'.format(h, t))
            break
     
    return  t,x,x_bind,x_nucl,sum_h



# @profile
def inf_model(t,x,x_bind,x_nucl,
        r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, 
        mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co,
        f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2,
        tau_death_Co, k_int_I1, k_int_I2, k_int_Co1,
        eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1,
        k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, 
        f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, 
        f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,
        k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo,
        endCo, startB1,startB2, n_bins_inf,age_bins,
        ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys):
    
    dxdt = np.zeros(len(x))  # Pre-allocate

    b = r_bleed * D
    
    # Species
    T = x[0]
    V1 = x[1]
    V2 = x[2]
    I1 = x[startI1:endI1+1]
    I2 = x[startI2:endI2+1]
    Co = x[startCo:endCo+1]
    NV = x[endCo+1]
    S = x[-1]
    
    # Pre-allocate all balances
    dxbind_dt = np.zeros(len(x_bind))
    dxnucl_dt = np.zeros(len(x_nucl))
    
    totI1 = np.sum(I1)
    totI2 = np.sum(I2)
    totCo = np.sum(Co)
    
    # Pre-allocate release scalars
    rel_v1 = 0
    rel_v2 = 0
    
    # Uninfected cells: infection and growth
    bind1 = k_bind_I1[0] * T * V1  # V1 infection rate counter
    bind2 = k_bind_I2[0] * T * V2  # V2 infection rate counter
    
    growth_lim = max(0, S / (K_s + S))
    
    growthT = mu_T * T * growth_lim
    dxdt[0] = growthT - k_deathT * T - (bind1 + bind2) + Cin * D / scaling - b * x[0]  # Target cells balance
    dxdt[endCo+1] += k_deathT * T
    
    dxdt[startI1] += bind1  # Infected cells balance
    dxdt[startI2] += bind2  # Infected cells balance
    
    dxbind_dt[startB1] += bind1  # Bound virus added to bound virus balance of infected cells of first age bin
    dxbind_dt[startB2] += bind2  # Bound virus added to bound virus balance of infected cells of first age bin
    
    # I1, I2: infection, death, and intracellular species balances
    j = np.arange(n_bins)
    age = age_bins


    # Growth
    growthI1 = mu_I1 * I1 * growth_lim
    growthI2 = mu_I2 * I2 * growth_lim
    dxdt[startI1 + j] += growthI1
    dxdt[startI2 + j] += growthI2

    # Death dependency on the number of DNA copies
    N1_pc=x_nucl[j] / (I1 + 1e-50)
    N2_pc=x_nucl[startB2 + j] / (I2 + 1e-50)
    
    N1_pc[N1_pc<np.e]=1
    N2_pc[N2_pc<np.e]=1
    k_deathI1 = k_deathT * (0.5 + 0.5 * np.tanh((tau_death_I1 - age) / 0.3)) + \
                k_death_I1 * np.log(N1_pc) * \
                (0.5 + 0.5 * np.tanh((age - tau_death_I1) / 0.3))

    k_deathI2 = k_deathT * (0.5 + 0.5 * np.tanh((tau_death_I2 - age) / 0.3)) + \
                k_death_I2 * np.log(N1_pc) * \
                    (0.5 + 0.5 * np.tanh((age - tau_death_I2) / 0.3))

    # I1 and I2 death
    deathI1 = k_deathI1 * I1
    deathI2 = k_deathI2 * I2
    dxdt[startI1 + j] -= deathI1
    dxdt[startI2 + j] -= deathI2
    dxdt[endCo+1] += deathI1.sum() + deathI2.sum()

    # V1, V2 bound to I1, I2: consumption for cell death and endocytosis
    Ein1 = k_int_I1 * x_bind[j]
    Ein2 = k_int_I2 * x_bind[startB2 + j]
    dxbind_dt[j] -= Ein1 + k_deathI1 * x_bind[j]
    dxbind_dt[startB2 + j] -= Ein2 + k_deathI2 * x_bind[startB2 + j ]
    dxnucl_dt[j] += Ein1 * eta_I1 - k_deathI1 * x_nucl[j] - k_d_N1 * x_nucl[j]
    dxnucl_dt[startB2 + j] += Ein2 * eta_I2 - k_deathI2 * x_nucl[startB2 + \
                 j] - k_d_N2 * x_nucl[startB2 + j]

    # Nuclear reactions in I1 and I2: BV DNA replication
    dxnucl_dt[j] += k_repl_I1 * x_nucl[j] * f_repl_I1[j]
    dxnucl_dt[startB2 + j] += k_repl_I2 * x_nucl[startB2 + j] * f_repl_I2[j]

    # Release from I1
    rel_v1 += np.sum(k_rel_I1 * I1 * f_rel_I1)
    rel_v2 += np.sum(k_rel_I1 * I1 * f_rel_I1 * randV2_I1)
    
    # Release from I2
    rel_v1 += np.sum(k_rel_I2 * I2 * f_rel_I2 * randV1_I2)
    rel_v2 += np.sum(k_rel_I2 * I2 * f_rel_I2)

    # Bleeding
    dxdt[startI1 + j] -= b * I1
    dxdt[startI2 + j] -= b * I2
    dxbind_dt[startB1 + j] -= b * x_bind[startB1 + j]
    dxbind_dt[startB2 + j] -= b * x_bind[startB2 + j]
    dxnucl_dt[startB1 + j] -= b * x_nucl[startB1 + j]
    dxnucl_dt[startB2 + j] -= b * x_nucl[startB2 + j]

    # Viral genome to nonviable cells
    dxnucl_dt[-2] += np.sum(k_deathI1 * x_nucl[j])
    dxnucl_dt[-1] += np.sum(k_deathI2 * x_nucl[startB2 + j])
    
    j = np.arange(n_bins_inf)
    # Co-infection: virus of a different type enters infected cell
    kv1 = k_bind_I2[j] * V1
    kv2 = k_bind_I1[j] * V2
    B1_new = I2[j] * kv1  # New virus 1 binding I2
    B2_new = I1[j] * kv2  # New virus 2 binding I1
    bind1 += B1_new.sum()  # Update counters of bound virus
    bind2 += B2_new.sum()  # Update counters of bound virus

    dxdt[startI1 + j] -= B2_new  # Remove new co-infected cells from I1 balance
    dxdt[startI2 + j] -= B1_new  # Remove new co-infected cells from I2 balance

    dxdt[startCo + ind[n_bins * j]] += B2_new  # Add new co-infected cells to Co balance
    dxdt[startCo + j] += B1_new  # Add new co-infected cells to Co balance

    v1_toCoB = kv2 * x_bind[j]  # V1 bound to I1 moving to co-infected cells
    v2_toCoB = kv1 * x_bind[startB2 + j]  # V2 bound to I2 moving to co-infected cells
    v1_toCoN = kv2 * x_nucl[j]  # E1 in I1 moving to co-infected cells
    v2_toCoN = kv1 * x_nucl[startB2 + j]  # E2 in I2 moving to co-infected cells

    dxbind_dt[j] -= v1_toCoB  # Move V1 bound to I1 that just got co-infected from B/I1 to B/Co balance
    dxbind_dt[startB2 + j] -= v2_toCoB  # Move V2 bound to I2 that just got co-infected from B/I2 to B/Co balance
    dxnucl_dt[j] -= v1_toCoN  # Move V1 bound to I1 that just got co-infected from B/I1 to B/Co balance
    dxnucl_dt[startB2 + j] -= v2_toCoN  # Move V2 bound to I2 that just got co-infected from B/I2 to B/Co balance

    dxbind_dt[startB1Co + j] += B1_new  # Update V1 bound to new co-infected cells
    dxbind_dt[ind_scB1Co[n_bins * j]] += v1_toCoB  # Update V1 bound to new co-infected cells
    dxnucl_dt[ind_scB1Co[n_bins * j]] += v1_toCoN  # Update E1 to new co-infected cells

    dxbind_dt[startB2Co + j] += v2_toCoB  # Update V2 bound to new co-infected cells
    dxnucl_dt[startB2Co + j] += v2_toCoN  # Update E2 bound to new co-infected cells
    dxbind_dt[ind_scB2Co[n_bins * j]] += B2_new  # Update V2 bound to new co-infected cells

    # Virus of the same type enters I1 and I2
    I1_uptake = k_bind_I1[j] * I1[j] * V1  # Calculate virus uptake
    I2_uptake = k_bind_I2[j] * I2[j] * V2
    dxbind_dt[j] += I1_uptake  # Add uptake to bound virus balance
    dxbind_dt[startB2 + j] += I2_uptake
    bind1 += I1_uptake.sum()  # Update counters
    bind2 += I2_uptake.sum()
        
        
    
    # Growth
    growthCo = mu_Co * Co * growth_lim
    dxdt[startCo:endCo+1] += growthCo

    rel_v1 = rel_v1 + np.sum(f_rel_Co * k_rel_Co * Co * (x_nucl[startB1Co:(startB1Co+len(Co))]) / \
            (x_nucl[startB1Co:(startB1Co+len(Co))] + x_nucl[startB2Co:(startB2Co+len(Co))] + 1e-30))
    rel_v2 = rel_v2 + np.sum(f_rel_Co * k_rel_Co * Co * (x_nucl[startB2Co:(startB2Co+len(Co))]) / \
            (x_nucl[startB1Co:(startB1Co+len(Co))] + x_nucl[startB2Co:(startB2Co+len(Co))] + 1e-30))

    N_Co_pc=(f_death_Co1 * x_nucl[startB1Co:(startB1Co+len(Co))] + \
             f_death_Co2 * x_nucl[startB2Co:(startB2Co+len(Co))]) / (Co + 1e-30)
    N_Co_pc[N_Co_pc<np.e]=1

    k_deathCo = k_deathT * (0.5 + 0.5 * np.tanh((tau_death_Co - age_co) / 0.3)) + \
            k_death_Co * np.log(N_Co_pc) * \
            (0.5 + 0.5 * np.tanh((age_co - tau_death_Co) / 0.3))

    
    deathCo = k_deathCo * Co
    dxdt[startCo:endCo+1] -= deathCo
    dxdt[endCo + 1] += deathCo.sum()

    # V1 binding to co-infected cells
    Co_uptakeV1 = k_bind_Co1 * Co * V1
    # V2 binding to co-infected cells
    Co_uptakeV2 = k_bind_Co2 * Co * V2

    Ein1 = k_int_Co1 * x_bind[startB1Co:(startB1Co+len(Co))]
    Ein2 = k_int_Co2 * x_bind[startB2Co:(startB2Co+len(Co))]

    dxbind_dt[startB1Co:(startB1Co+len(Co))] += Co_uptakeV1 - k_deathCo * \
        x_bind[startB1Co :(startB1Co+len(Co))] - Ein1
    dxbind_dt[startB2Co:(startB2Co+len(Co))] += Co_uptakeV2 - k_deathCo * \
        x_bind[startB2Co:(startB2Co+len(Co))] - Ein2

    dxnucl_dt[startB1Co:(startB1Co+len(Co))] += Ein1 * eta_Co1 - \
        k_deathCo * x_nucl[startB1Co:(startB1Co+len(Co))] - k_d_N1 *\
            x_nucl[startB1Co:(startB1Co+len(Co))]
    dxnucl_dt[startB2Co:(startB2Co+len(Co))] += Ein2 * eta_Co2 - \
        k_deathCo * x_nucl[startB2Co:(startB2Co+len(Co))] - k_d_N2 * \
            x_nucl[startB2Co:(startB2Co+len(Co))]

    bind1 = bind1 + Co_uptakeV1.sum()
    bind2 = bind2 + Co_uptakeV2.sum()

    # Nuclear reactions in co-infected cells: BV DNA replication
    dxnucl_dt[startB1Co:(startB1Co+len(Co))] += k_repl_Co1 * \
        x_nucl[startB1Co:(startB1Co+len(Co))] * f_repl_Co1
    dxnucl_dt[startB2Co:(startB2Co+len(Co))] += k_repl_Co2 * \
        x_nucl[startB2Co:(startB2Co+len(Co))] * f_repl_Co2

    # Bleeding
    dxdt[startCo:(startCo+len(Co))] -= b * Co
    dxbind_dt[startB1Co:(startB1Co+len(Co))] -= b * x_bind[startB1Co:(startB1Co+len(Co))]
    dxnucl_dt[startB1Co:(startB1Co+len(Co))] -= b * x_nucl[startB1Co:(startB1Co+len(Co))]
    dxbind_dt[startB2Co:(startB2Co+len(Co))] -= b * x_bind[startB2Co:(startB2Co+len(Co))]
    dxnucl_dt[startB2Co:(startB2Co+len(Co))] -= b * x_nucl[startB2Co:(startB2Co+len(Co))]

    # Viral genome to nonviable cells
    dxnucl_dt[-2] += np.sum(k_deathCo * x_nucl[startB1Co:(startB1Co+len(Co))])
    dxnucl_dt[-1] += np.sum(k_deathCo * x_nucl[startB2Co:(startB2Co+len(Co))])
    
        
    # Virions balance
    dxdt[1] = -bind1 - k_degrBV_1 * V1 + rel_v1 - D * x[1]  # BV ITR/GOI (#/mL)
    dxdt[2] = -bind2 - k_degrBV_2 * V2 + rel_v2 - D * x[2]  # BV rep/cap (#/mL)
    
    # Nonviable cells: outlet, lysis, viral genome degradation
    dxdt[endCo+1] -= (k_lys + b) * NV
    dxnucl_dt[-2] -= (k_lys + b + k_d_N1) * x_nucl[-2]
    dxnucl_dt[-1] -= (k_lys + b + k_d_N2) * x_nucl[-1]
    
    # Substrate consumption
    dxdt[-1] = D * (Sin - S) - (Ys_T * T + Ys_I1 * totI1 + Ys_I2 * totI2 + Ys_C * totCo) * S / (1e-2 / scaling + S)
    
    return dxdt,dxbind_dt,dxnucl_dt