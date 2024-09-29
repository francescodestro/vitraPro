"""
vitraPro: Simulation of viral transduction and propagation in suspension cultures

    v1.0 July 24, 2024. F. Destro
    
    Reference:  
        Destro F. and R.D. Braatz (2024). Efficient simulation of viral
        transduction and propagation for biomanufacturing. ACS
        Synthetic Biology.

Outputs:

    Scalar states:
        T:              Uninfected cell concentration time profile - [cell/mL]
        V1:             Virion 1 concentration time profile - [PFU/mL]
        NV:             Nonviable cell concentration time profile  - [cell/mL]
        S:              Substrate concentration time profile - [nmol/mL]
        N1_pc_avg:      Time profile of avg concentration of virus 1 genome in 
                        infected cells - [vg/cell] (per cell basis)

    States distributed with respect to one infection age:
        I1:             Time profile of concentration of cells infected by 
                        virus 1 - [cell/mL] 
        B1:             Time profile of concentration of virus 1 bound to infected cells 
                        - [vg/mL] (total conc. in system)
        B1_pc:          Time profile of concentration of virus 1 bound to infected cells 
                        - [vg/mL] (total conc. in system)
        N1:             Time profile of concentration of virus 1 genome in 
                        nucleus of infected cells - [vg/mL] (total conc. in system)
        N1_pc:          Time profile of concentration of virus 1 genome in 
                        nucleus of infected cells - [vg/mL] (total conc. in system)

"""

import numpy as np
from def_parameters import parameters
from main_functions import main_RK23
from sample_figures import plot_sample_figures
import warnings

########################################################################
#   
#   Simulation inputs setup
#
########################################################################

t_final = 24 * 10  # simulation duration [h]
C0_v = 1.5e6  # viable cell density at t=0 [cell/mL]
C0_nv = 0  # nonviable cell density at t=0 [cell/mL]
S0 = 15 * 1e3  # substrate concentration [nmol/mL]

# Inoculation of virions
MOI_1 = 1  # [PFU/viable cell]

# Inoculation of infected cells
I0_1 = 0 * C0_v / 100  # [#/mL] - concentration of inoculated cells infected by virus 1
N_I0_1 = 5e4 * I0_1  # [vg/cell] - viral genome copy number in inoculated cells infected by virus 1
age_I0_1 = 40  # [hpi] - infection age of inoculated cells infected by virus 

# Discretization parameters
Dtau = 0.2  # [hpi] mesh spacing
age_viab = 120  # [hpi] - maximum infection age in the mesh
age_end_inf = 10  # [hpi] - infection age beyond which no viral binding occurs

# Inlet conditions (see Table 1 in Destro and Braatz, 2024)
Cin = 3e6  # [cell/mL] - viable uninfected cells in feed
Din = 0.01  # [1/h] - dilution rate
Sin = 15 * 1e3  # [nmol/mL] - substrate concentration in feed
r_bleed = 1  # bleeding ratio (0 ≤ r < 1 for perfusion, r = 1 for continuous)

# Inlet profiles
Dt = 1  # [h] control and sampling interval
tt = np.arange(0, t_final + Dt, Dt)
CCin = np.ones(len(tt)) * Cin
DD = np.ones(len(tt)) * Din
r_bleed_vect = np.ones(len(tt)) * r_bleed
Sin_vect = np.ones(len(tt)) * Sin

# Numerical scheme
numeric_scheme = 'RK23'

generate_figures = True

# Simulation
scaling = 5e7
p = parameters(age_viab, age_end_inf, Dtau, scaling)  

n_bins = p['n_bins']
n_bins_inf = p['n_bins_inf']
t = 0
V0_1 = C0_v * MOI_1  # PFU/mL
x0 = np.array([C0_v, V0_1] + [0] * n_bins + [C0_nv, S0]) / scaling
x0_bind = np.zeros(n_bins)
x0_nucl = np.hstack((x0_bind, 0))
sum_h = 0

# Preallocate solution vectors
xx = np.zeros((len(tt), len(x0)))
xx_bind = np.zeros((len(tt), n_bins))
xx_nucl = np.zeros((len(tt), n_bins + 1))
xx[0, :] = x0 * scaling

Sin_vect = Sin_vect / scaling

# Inoculation of infected cells
bin_inoculation_1 = round(age_I0_1 / Dtau)
x0[p['startI1'] + bin_inoculation_1] = I0_1 / scaling
x0_nucl[p['startB1'] + bin_inoculation_1] = N_I0_1 / scaling

for i in range(1, len(tt)):

    D = DD[i]
    Cin = CCin[i]
    r_bleed = r_bleed_vect[i]
    Sin = Sin_vect[i]

    t_out, x, x_bind, x_nucl, sum_h = main_RK23([t, t + Dt], x0, x0_bind, x0_nucl, 
                                           D, r_bleed, Cin, Sin, p, sum_h, scaling)

    # Update for the next step
    x0 = x
    x0_bind = x_bind
    x0_nucl = x_nucl
    t = t_out

    # Save new results
    index = i
    xx[index, :] = x * scaling
    xx_bind[index, :] = x_bind * scaling
    xx_nucl[index, :] = x_nucl * scaling


# Simulation outputs
T = xx[:, 0]  # [cell/mL] - scalar
V1 = xx[:, 1]  # [PFU/mL] - scalar
NV = xx[:, -2]  # [cell/mL] - scalar
S = xx[:, -1]  # [nmol/mL] - scalar

warnings.filterwarnings("ignore")

I1 = xx[:, p['startI1']:p['endI1']+1]  # [cell/mL] - distribution
B1 = xx_bind[:, p['startB1']:p['endB1']+1]  # [vg/mL] - distribution
B1_pc = B1 / I1  # [vg/cell] - distribution
N1 = xx_nucl[:, p['startB1']:p['endB1']+1]  # [vg/mL] - distribution
N1_pc = N1 / I1  # [vg/cell] - distribution

N1_pc_avg = np.sum(N1, axis=1) / np.sum(I1, axis=1)  # [vg/cell] - scalar

# Plots
if generate_figures:
    plot_sample_figures(tt, T, V1, I1, NV, S, N1_pc_avg, p, N1_pc)

