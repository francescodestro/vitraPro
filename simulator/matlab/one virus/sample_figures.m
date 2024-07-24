% Generate sample figures
col='b'; % line color
fs=18; % fontsize

% Uninfected cells concentration
figure(1)
hold on 
plot(tt/24,T,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Uninfected cells [cell/mL]')
xlabel('Time [d]')
box on

% Virion 1 concentration
figure(2)
hold on
plot(tt/24,V1,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
xlabel('Time [d]')
ylabel('Virion 1 [PFU/mL]')
box on

% Cells infected by virus 1: total concentration
figure(4)
hold on
plot(tt/24,sum(I1,2),'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs)
ylabel('Cells infected by virus 1 [cell/mL]')
box on

% Nonviable cells concentration
figure(7)
hold on
plot(tt/24,NV,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Nonviable cells [cell/mL]')
box on

% Substrate concentration
figure(8)
hold on
plot(tt/24,S,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Substrate [nmol/mL]')
box on

% Average virus 1 copy number in nucleus of infected cells
figure(9)
hold on 
plot(tt/24,N1_pc_avg,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Avg virus 1 in nucleus of I1 [vg/cell]')
xlabel('Time [d]')
box on

%% Examples of plots of distributions
[~,t_plot]=min(abs(tt-24)); % tt(25) = 24 h with default sampling interval

% Infection age distribution of I1 at t = 24 h
figure(13)
hold on 
plot(p.age_bins,I1(t_plot,:),'o','linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Age distribution of I1 [cell/mL]')
xlabel('Infection age [hpi]')
box on

% Virus 1 copy number in nucleus of infected cells at t = 24 h
figure(15)
hold on 
plot(p.age_bins,N1_pc(t_plot,:),'o','linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Distribution of virus 1 in nucleus of I1 [vg/cell]')
xlabel('Infection age [hpi]')
box on
