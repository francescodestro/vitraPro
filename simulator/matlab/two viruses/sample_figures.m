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

% Virion 2 concentration
figure(3)
hold on
plot(tt/24,V2,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
xlabel('Time [d]')
ylabel('Virion 2 [PFU/mL]')
box on

% Cells infected by virus 1: total concentration
figure(4)
hold on
plot(tt/24,sum(I1,2),'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs)
ylabel('Cells infected by virus 1 [cell/mL]')
box on

% Cells infected by virus 2: total concentration
figure(5)
hold on
plot(tt/24,sum(I2,2),'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs)
ylabel('Cells infected by virus 2 [cell/mL]')
box on

% Coinfected cells: total concentration
figure(6)
hold on
plot(tt/24,sum(Co,2),'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Coinfected cells [cell/mL]')
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

% Average virus 2 copy number in nucleus of infected cells
figure(10)
hold on 
plot(tt/24,N2_pc_avg,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Avg virus 2 in nucleus of I2 [vg/cell]')
xlabel('Time [d]')
box on

% Average virus 1 copy number in nucleus of infected cells
figure(11)
hold on 
plot(tt/24,N1_Co_pc_avg,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Avg virus 1 in nucleus of coinfected [vg/cell]')
xlabel('Time [d]')
box on

% Average virus 2 copy number in nucleus of infected cells
figure(12)
hold on 
plot(tt/24,N2_Co_pc_avg,'linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Avg virus 2 in nucleus of coinfected [vg/cell]')
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

% Infection age distribution of I2 at t = 24 h
figure(14)
hold on 
plot(p.age_bins,I2(t_plot,:),'o','linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Age distribution of I2 [cell/mL]')
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

% Virus 2 copy number in nucleus of infected cells at t = 24 h
figure(16)
hold on 
plot(p.age_bins,N2_pc(t_plot,:),'o','linewidth',1.5,'Color',col)
set(gca,'linewidth',1.5,'fontsize',fs) 
ylabel('Distribution of virus 2 in nucelus of I2 [vg/cell]')
xlabel('Infection age [hpi]')
box on

% Infection age distribution of coinfected cells at t = 24 h
figure(17)
max_age=60; % [hpi] - max age for x and y axis in 3D plot
age_step=20; % [hpi] - tick for x and y axis in 3D plot

co=Co(t_plot,:);
for i1=1:n_bins
    for i2=1:n_bins
        index=p.ind(i2+n_bins*(i1-1));
        if index >0 
            co_plot(i1,i2)=co(index);
        end
    end
end
h=bar3(co_plot);
view([-90 30 60])
for k=1:length(h)
    h(k).LineStyle='none';
end
set(gca,'GridAlpha',0.1)  
axis_tick=p.age_bins;
step_size=age_step/p.Dtau+1;
set(gca,'fontsize',28,'LineWidth',2,'xtick',...
    0:step_size-1:length(co),'xticklabel',...
    axis_tick(1):axis_tick(step_size):(axis_tick(end)+p.Dtau),...
    'ytick',0:step_size-1:length(co),'yticklabel',...
    axis_tick(1):axis_tick(step_size):(axis_tick(end)+p.Dtau),...
    'xlim',[0 max_age/p.Dtau],'ylim',[0 max_age/p.Dtau])
ylabel('Infection age 1 [hpi]')
xlabel('Infection age 2 [hpi]')
zlabel('Coinfected cells [cell/mL]')

% Viral genome distribution in nucleus of coinfected cells at t = 24 h
figure(18)
max_age=60; % [hpi] - max age for x and y axis in 3D plot
age_step=20; % [hpi] - tick for x and y axis in 3D plot

co=Co(t_plot,:);
n1_co=xx_nucl(t_plot,p.startB1Co:p.endB1Co);

for i1=1:n_bins
    for i2=1:n_bins
        index=p.ind(i2+n_bins*(i1-1));
        if index >0 
            co_plot(i1,i2)=co(index);
            n1_co_plot(i1,i2)=n1_co(index);
        end
    end
end
h=bar3((n1_co_plot)./co_plot);
view([-90 30 60])
for k=1:length(h)
    h(k).LineStyle='none';
end
set(gca,'GridAlpha',0.1)  
axis_tick=p.age_bins;
step_size=age_step/p.Dtau+1;
set(gca,'fontsize',28,'LineWidth',2,'xtick',...
    0:step_size-1:length(co),'xticklabel',...
    axis_tick(1):axis_tick(step_size):(axis_tick(end)+p.Dtau),...
    'ytick',0:step_size-1:length(co),'yticklabel',...
    axis_tick(1):axis_tick(step_size):(axis_tick(end)+p.Dtau),...
    'xlim',[0 max_age/p.Dtau],'ylim',[0 max_age/p.Dtau])
ylabel('Infection age 1 [hpi]')
xlabel('Infection age 2 [hpi]')
zlabel('Viral genome 1 [vg/cell]')


