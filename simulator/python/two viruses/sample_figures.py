import matplotlib.pyplot as plt
import numpy as np

def plot_sample_figures(tt,T,V1,V2,I1,I2,Co,NV,S,N1_pc_avg,N2_pc_avg,N1_Co_pc_avg,
                        N2_Co_pc_avg,N1_pc,N2_pc,n1_co,p):

    # Define line color and font size
    col = 'b'  # line color
    fs = 18    # fontsize
    

    # Uninfected cells concentration
    plt.figure(1)
    plt.plot(tt/24, T, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Uninfected cells [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Virion 1 concentration
    plt.figure(2)
    plt.plot(tt/24, V1, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Virion 1 [PFU/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Virion 2 concentration
    plt.figure(3)
    plt.plot(tt/24, V2, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Virion 2 [PFU/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Cells infected by virus 1: total concentration
    plt.figure(4)
    plt.plot(tt/24, np.sum(I1,1), linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Cells infected by virus 1 [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)

    # Cells infected by virus 2: total concentration
    plt.figure(5)
    plt.plot(tt/24, np.sum(I2,1), linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Cells infected by virus 2 [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Coinfected cells: total concentration
    plt.figure(6)
    plt.plot(tt/24, np.sum(Co, 1), linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Coinfected cells [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Nonviable cells concentration
    plt.figure(7)
    plt.plot(tt/24, NV, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Nonviable cells [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Substrate concentration
    plt.figure(8)
    plt.plot(tt/24, S, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Substrate [nmol/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Average virus 1 copy number in nucleus of infected cells
    plt.figure(9)
    plt.plot(tt/24, N1_pc_avg, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Avg virus 1 in nucleus of I1 [vg/cell]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)

    # Average virus 2 copy number in nucleus of infected cells
    plt.figure(10)
    plt.plot(tt/24, N2_pc_avg, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Avg virus 2 in nucleus of I2 [vg/cell]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Average virus 1 copy number in nucleus of coinfected cells
    plt.figure(11)
    plt.plot(tt/24, N1_Co_pc_avg, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Avg virus 1 in nucleus of coinfected [vg/cell]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Average virus 2 copy number in nucleus of coinfected cells
    plt.figure(12)
    plt.plot(tt/24, N2_Co_pc_avg, linewidth=1.5, color=col)
    plt.xlabel('Time [d]')
    plt.ylabel('Avg virus 2 in nucleus of coinfected [vg/cell]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
                                                
    t_plot = 24
    
    # Infection age distribution of I1 at t = 24 h
    plt.figure(13)
    plt.plot(p['age_bins'], I1[t_plot, :], 'o', linewidth=1.5, color=col)
    plt.xlabel('Infection age [hpi]')
    plt.ylabel('Age distribution of I1 [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Infection age distribution of I2 at t = 24 h
    plt.figure(14)
    plt.plot(p['age_bins'], I2[t_plot, :], 'o', linewidth=1.5, color=col)
    plt.xlabel('Infection age [hpi]')
    plt.ylabel('Age distribution of I2 [cell/mL]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Virus 1 copy number in nucleus of infected cells at t = 24 h
    plt.figure(15)
    plt.plot(p['age_bins'], N1_pc[t_plot, :], 'o', linewidth=1.5, color=col)
    plt.xlabel('Infection age [hpi]')
    plt.ylabel('Distribution of virus 1 in nucleus of I1 [vg/cell]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # Virus 2 copy number in nucleus of infected cells at t = 24 h
    plt.figure(16)
    plt.plot(p['age_bins'], N2_pc[t_plot, :], 'o', linewidth=1.5, color=col)
    plt.xlabel('Infection age [hpi]')
    plt.ylabel('Distribution of virus 2 in nucleus of I2 [vg/cell]')
    plt.gca().tick_params(width=1.5, labelsize=fs)
    plt.box(True)
    
    # 3d plots: define parameters
    t_plot = 24  # tt(25) = 24 h with default sampling interval
    max_age = 60  # [hpi] - max age for x and y axis 
    age_step = 20  # [hpi] - tick for x and y axis
    
    # Create 3D bar plot for infection age distribution
    co = Co[t_plot, :]
    n_bins=p['n_bins']
    co_plot = np.zeros((n_bins, n_bins))
    for i1 in range(n_bins):
        for i2 in range(n_bins):
            index = p['ind'][i2 + n_bins * (i1 - 1)]
            if index > 0:
                co_plot[i1, i2] = co[index]
    # fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111, projection='3d')
    # xpos, ypos = np.meshgrid(range(n_bins), range(n_bins))
    # xpos = xpos.flatten()
    # ypos = ypos.flatten()
    # zpos = np.zeros_like(xpos)
    # dx = dy = 0.8 * np.ones_like(zpos)
    # dz = co_plot.flatten()
    
    # bars = ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b')
    
    # # ax.view_init(-90, 30)
    # ax.grid(alpha=0.1)
    # axis_tick = p['age_bins']
    # step_size = int(age_step / p['Dtau']) + 1
    # ax.set_xlabel('Infection age 1 [hpi]')
    # ax.set_ylabel('Infection age 2 [hpi]')
    # ax.set_zlabel('Co')
    # ax.set_xticks(np.arange(0, len(co), step_size))
    # ax.set_xticklabels(axis_tick[::step_size])
    # ax.set_yticks(np.arange(0, len(co), step_size))
    # ax.set_yticklabels(axis_tick[::step_size])
    # ax.set_xlim(0, max_age / p['Dtau'])
    # ax.set_ylim(0, max_age / p['Dtau'])

    # # reate 3D bar plot for viral genome distribution in nucleus of coinfected cells at t = 24 h
    # # n1_co = xx_nucl[t_plot, p['startB1Co']:p['endB1Co']]
    # n1_co_plot = np.zeros((n_bins, n_bins))
    # for i1 in range(n_bins):
    #     for i2 in range(n_bins):
    #         index = p['ind'][i2 + n_bins * (i1 - 1)]
    #         if index > 0:
    #             n1_co_plot[i1, i2] = n1_co[index]
    
    # plt.figure(18)
    # h = plt.bar3d(range(n_bins), range(n_bins), np.zeros_like(n1_co_plot).ravel(), 1, 1, (n1_co_plot / co_plot).ravel())
    # plt.gca().view_init(elev=-90, azim=30)
    # plt.gca().grid(alpha=0.1)
    # plt.xticks(np.arange(0, len(co), step_size), axis_tick[::step_size])
    # plt.yticks(np.arange(0, len(co), step_size), axis_tick[::step_size])
    # plt.xlim(0, max_age / p['Dtau'])
    # plt.ylim(0, max_age / p['Dtau'])
    # plt.gca().tick_params(axis='both', which='major', labelsize=28, width=2)
    # plt.xlabel('Infection age 1 [hpi]')
    # plt.ylabel('Infection age 2 [hpi]')
    # plt.zlabel('Viral genome 1 [vg/cell]')
    
    plt.show()




