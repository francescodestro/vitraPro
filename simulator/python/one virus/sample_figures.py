import matplotlib.pyplot as plt
import numpy as np

def plot_sample_figures(tt, T, V1, I1, NV, S, N1_pc_avg, p, N1_pc):
    col = 'b'  # line color
    fs = 18  # fontsize

    # Uninfected cells concentration
    plt.figure(1)
    plt.plot(tt / 24, T, linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Uninfected cells [cell/mL]', fontsize=fs)
    plt.xlabel('Time [d]', fontsize=fs)
    plt.box(True)

    # Virion 1 concentration
    plt.figure(2)
    plt.plot(tt / 24, V1, linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.xlabel('Time [d]', fontsize=fs)
    plt.ylabel('Virion 1 [PFU/mL]', fontsize=fs)
    plt.box(True)

    # Cells infected by virus 1: total concentration
    plt.figure(4)
    plt.plot(tt / 24, np.sum(I1, axis=1), linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Cells infected by virus 1 [cell/mL]', fontsize=fs)
    plt.box(True)

    # Nonviable cells concentration
    plt.figure(7)
    plt.plot(tt / 24, NV, linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Nonviable cells [cell/mL]', fontsize=fs)
    plt.box(True)

    # Substrate concentration
    plt.figure(8)
    plt.plot(tt / 24, S, linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Substrate [nmol/mL]', fontsize=fs)
    plt.box(True)

    # Average virus 1 copy number in nucleus of infected cells
    plt.figure(9)
    plt.plot(tt / 24, N1_pc_avg, linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Avg virus 1 in nucleus of I1 [vg/cell]', fontsize=fs)
    plt.xlabel('Time [d]', fontsize=fs)
    plt.box(True)

    # Find closest time point to t = 24 h
    t_plot_idx = np.argmin(np.abs(tt - 24))

    # Infection age distribution of I1 at t = 24 h
    plt.figure(13)
    plt.plot(p['age_bins'], I1[t_plot_idx, :], 'o', linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Age distribution of I1 [cell/mL]', fontsize=fs)
    plt.xlabel('Infection age [hpi]', fontsize=fs)
    plt.box(True)

    # Virus 1 copy number in nucleus of infected cells at t = 24 h
    plt.figure(15)
    plt.plot(p['age_bins'], N1_pc[t_plot_idx, :], 'o', linewidth=1.5, color=col)
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.ylabel('Distribution of virus 1 in nucleus of I1 [vg/cell]', fontsize=fs)
    plt.xlabel('Infection age [hpi]', fontsize=fs)
    plt.box(True)

    # Display all figures
    plt.show()

