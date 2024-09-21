% Simulation of viral transduction and propagation in suspension cultures
%   This script generates Figures 6 and 7
clear
clc

%% Case study 3: Initialization and simulation

% Initialization
cols=[{[0 169 192]/255} {[51 74 159]/255} {[173 76 157]/255}  {[0.7 0.7 0.7]} {'k'} {[0.4660, 0.6740, 0.1880]}];
lines=[{':'},{'-'},{'-.'}];
lws=[3.5 3 3.5];
ylims_I=[5e3 5e3 5e3 .5e3]*10;
ylims_DNA=[13e7 13e7 13e7 500];
lw=lws(1);
line=lines{2};
col=[0, 0.4470, 0.7410];
scaling=5e7;
pos=[2079 519 473 245];

% numeric_scheme='RK45';
numeric_scheme='RK23';

% Simulation parameters
t_final=24*10; % simulation duration [h]
C0_v=1.5e6; % viable cell density at t=0 [#/mL]
C0_nv=0; % nonviable cell density at t=0 [#/mL]
S0=15*1e3; % substrate concentration [nmol/mL]

% Inoculation with infected cells
I0_1=C0_v/100;
N_I0_1=5e4*I0_1;
age_I0_1=40; %hpi
I0_2=C0_v/100/2;
N_I0_2=5e4*I0_2;
age_I0_2=40; %hpi

% Inoculation with virions
MOI_1=0;
MOI_2=0;

% Inlet parameters
Cin=3e6;
Din=0.01;
Sin=15*1e3/scaling;
r_bleed=1;

% Discretization parameters
Dtau=.2; % mesh grid
age_end_inf=20;
age_viab=140;
Dt = 1; % h: control interval

p=def_parameters(age_viab,age_end_inf,Dtau,Dt,C0_v,MOI_1,MOI_2,scaling);

n_bins=p.n_bins;
n_bins_inf=p.n_bins_inf;
t=0;
V0_1=C0_v*MOI_1; % #/mL
V0_2=C0_v*MOI_2; % #/mL
x0=[C0_v V0_1 V0_2 zeros(1,n_bins*2+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1)) C0_nv S0]/scaling;

x0_bind=zeros(1,p.endB2Co);
x0_nucl=x0_bind;

sum_h=0;
    
% Preallocate solution vectors
tt=0:Dt:t_final;
xx=zeros(length(tt),length(x0));
xx_bind=zeros(length(tt),p.endB2Co);
xx_nucl=xx_bind;
xx(1,:)=x0*scaling;
    
CCin=ones(1,length(tt))*Cin;
DD=ones(1,length(tt))*Din;
r_bleed_vect=ones(1,length(tt))*r_bleed;

% Inoculation with infected cells
bin_inoculation_1=round(age_I0_1/Dtau);
x0(p.startI1+bin_inoculation_1)=I0_1/scaling;
x0_nucl(p.startB1+bin_inoculation_1)=N_I0_1/scaling;
bin_inoculation_2=round(age_I0_2/Dtau);
x0(p.startI2+bin_inoculation_2)=I0_2/scaling;
x0_nucl(p.startB2+bin_inoculation_2)=N_I0_2/scaling;

% Simulation
for i = 2:length(tt)

    D=DD(i);
    Cin=CCin(i);
    r_bleed=r_bleed_vect(i);

    switch numeric_scheme

        case 'RK23'
            [t_out,x,x_bind,x_nucl,sum_h] = main_RK23([t t+Dt], x0, x0_bind, ...
                x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling);

        case 'RK45'
            [t_out,x,x_bind,x_nucl,sum_h] = main_RK45([t t+Dt], x0, x0_bind, ...
                x0_nucl, D, r_bleed, Cin, Sin, p, sum_h, scaling);

    end

    % initialize following step
    x0=x(end,:);
    x0_bind=x_bind(end,:);
    x0_nucl=x_nucl(end,:);
    
    t=t_out(end);

    % save new results
    index=i;
    xx(index,:)=x(end,:)*scaling;
    xx_bind(index,:)=x_bind(end,:)*scaling;
    xx_nucl(index,:)=x_nucl(end,:)*scaling;

end

%% Figure 6
I1=xx(:,p.startI1:p.endI1);
I2=xx(:,p.startI2:p.endI2);
N1=xx_nucl(:,p.startB1:p.endB1);
Co=xx(:,p.startCo:p.endCo);

figure(1)
hold on
plot(tt/24,sum(I1,2),'linewidth',lws(1),'Color',cols{2},'LineStyle',lines{1})
plot(tt/24,sum(I2,2),'linewidth',lws(2),'Color',cols{3},'LineStyle',lines{2})
plot(tt/24,sum(Co,2),'linewidth',lws(3),'Color',cols{4},'LineStyle',lines{3})
legend('Virus 1','Virus 2','Coinfected')
set(gca,'linewidth',2,'fontsize',22,'ylim',[0 2e6]) %,'xticklabel',[])
set(gcf,'Position',pos)
xlabel('Time [d]')
ylabel('Infected [cell/mL]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(2)
hold on
plot(tt/24,xx(:,2),'linewidth',lws(1),'Color',cols{2},'LineStyle',lines{1})
plot(tt/24,xx(:,3),'linewidth',lws(2),'Color',cols{3},'LineStyle',lines{2})
set(gca,'linewidth',2,'fontsize',22,'ylim',[0 1e8],'xtick',0:2:10) %,'xticklabel',[])
set(gcf,'Position',pos)
ylabel('Virion [PFU/mL]')
xlabel('Time [d]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(4)
hold on 
plot(tt/24,xx(:,1),'linewidth',lw,'Color',cols{1},'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) %,'ylim',[0 4e6]) %,'xticklabel',[])
set(gcf,'Position',pos)
ylabel('Uninfected [cell/mL]')
xlabel('Time [d]')
box on

figure(5)
hold on 
plot(tt/24,sum(xx_nucl(:,p.startB1Co:p.endB1Co),2)./sum(Co,2),...
    'linewidth',lws(1),'Color',cols{2},'LineStyle',lines{1})
hold on
plot(tt/24,sum(xx_nucl(:,p.startB2Co:p.endB2Co),2)./sum(Co,2),...
    'linewidth',lws(2),'Color',cols{3},'LineStyle',lines{2})
set(gca,'linewidth',2,'fontsize',22) %,'ylim',[0 4e6]) %,'xticklabel',[])
set(gcf,'Position',[1819        -174         503         245])
ylabel([{'Avg viral genome in'},{'coinfected [vg/cell]'}])
xlabel('Time [d]')
box on

DNA1=xx_nucl(:,p.startB1:p.endB1);
DNA2=xx_nucl(:,p.startB2:p.endB2);

t_plot=[24 48 72 240]+1;
pos2=[ 323   454   560   150];

% infection age distribution (7E)
for i=1:length(t_plot)
    figure(5+i)
    hold on
    ylim_I=ylims_I(i);
    plot(p.age_bins,I1(t_plot(i),:),'linewidth',lws(1),'Color',cols{2},'LineStyle',lines{1})
    plot(p.age_bins,I2(t_plot(i),:),'linewidth',lws(2),'Color',cols{3},'LineStyle',lines{2})
    set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],'ylim',[0 ylim_I],'fontsize',22,'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
    box on
    set(gcf,'Position',pos2);
end
set(gca,'linewidth',2,'xtick',0:20:140,'xticklabel',0:20:140,'fontsize',22,'xlim',[0 100],'YScale','linear') %,'xticklabel',[])

% % DNA/mL
% for i=1:length(t_plot) 
%     figure(13+i)
%     hold on
%     plot(p.age_bins,DNA1(t_plot(i),:),'linewidth',lws(1),'Color',cols{1},'LineStyle',lines{1})
%     plot(p.age_bins,DNA2(t_plot(i),:),'linewidth',lws(2),'Color',cols{2},'LineStyle',lines{2})
%     set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],'fontsize',22,...
%         'xlim',[0 100],'YScale','linear','ylim',[0 ylims_DNA(i)])
%     box on
%     set(gcf,'Position',pos2);
%     
% end
% set(gca,'linewidth',2,'xtick',0:20:140,'xticklabel',0:20:140,'fontsize',22,'xlim',[0 100],'YScale','linear','ylim',[0 ylims_DNA(i)]) %'ylim',[0 2.5e8],,'xticklabel',[])

% Save figures for manuscript
% figure(1)
% set(gca,'xtick',0:2:10)
% saveas(gcf,'I_CS2.png')
% 
% figure(2)
% set(gca,'ylim',[0 8e7])
% set(gca,'xtick',0:2:10)
% saveas(gcf,'V_CS2.png')
% 
% figure(4)
% set(gca,'xtick',0:2:10)
% saveas(gcf,'T_CS2.png')
% 
% figure(5)
% set(gca,'xtick',0:2:10)
% saveas(gcf,'DNA_CS2ansf.png')
% 
% figure(6)
% set(gca,'ylim',[0 3e4])
% set(gca,'ytick',[0 2e4])
% text(55, 3e4*.8,'1 d post inoculation','FontSize',20)
% saveas(gcf,'age_CS2_1.png')
% 
% figure(7)
% set(gca,'ylim',[0 3e4])
% set(gca,'ytick',[0 2e4])
% text(55, 3e4*.8,'2 d post inoculation','FontSize',20)
% saveas(gcf,'age_CS2_2.png')
% 
% figure(8)
% set(gcf,'Position',[ 729   458   565   122])
% set(gca,'ylim',[0 3e4])
% set(gca,'ytick',[0 2e4])
% text(55, 3e4*.8,'3 d post inoculation','FontSize',20)
% saveas(gcf,'age_CS2_3.png')
% 
% figure(9)
% set(gcf,'Position',[  861   202   565   139])
% set(gca,'ylim',[0 4000])
% set(gca,'ytick',[0 2000],'yticklabel',[0 0.2])
% text(55, 4000*.8,'10 d post inoculation','FontSize',20)
% saveas(gcf,'age_CS2_4.png')


%% Figure 7: 3D plots

% Figure 7a
startB1Co=p.startB1Co;
endB1Co=p.endB1Co;
startB2Co=p.startB2Co;
endB2Co=p.endB2Co;
ind=p.ind;

for t_3Dplot=t_plot(2) 
    figure
    co=Co(t_3Dplot,:);
    b1_co=xx_nucl(t_3Dplot,startB1Co:endB1Co);
    b2_co=xx_nucl(t_3Dplot,startB2Co:endB2Co);
    co_plot=zeros(n_bins);
    b1_co_plot=co_plot;
    b2_co_plot=co_plot;
    
    for i1=1:n_bins
        for i2=1:n_bins
            index=ind(i2+n_bins*(i1-1));
            if index >0 
                co_plot(i1,i2)=co(index);
                b1_co_plot(i1,i2)=b1_co(index);
                b2_co_plot(i1,i2)=b2_co(index);
            end
        end
    end
    h=bar3(co_plot);
    view([-150 30 90])
    for k=1:length(h)
        h(k).LineStyle='none';
    end
    set(gca,'GridAlpha',0.1)
    if t_3Dplot ==t_plot(2)
        axis_tick=p.age_bins;
        step_size=51;
        set(gca,'fontsize',28,'LineWidth',2,'xtick',...
            0:step_size-1:length(co),'xticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'ytick',0:step_size-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'xlim',[0 251],'ylim',[0 251])
    else
        axis_tick=p.age_bins;
        step_size_y=101;
        step_size_x=201;
        set(gca,'fontsize',28,'LineWidth',2,'ytick',...
            0:step_size_y-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size_y):axis_tick(end)+Dtau,...
            'xtick',0:step_size_x-1:800,'xticklabel',...
            axis_tick(1):axis_tick(step_size_x):axis_tick(end)+Dtau,...
            'xlim',[0 601],'ylim',[0 601])
    end
    set(gcf,'Position',[348   259   798   586])
    colormap default
end

% Figure 7b
startB1Co=p.startB1Co;
endB1Co=p.endB1Co;
startB2Co=p.startB2Co;
endB2Co=p.endB2Co;
ind=p.ind;

for t_3Dplot=t_plot(2)
    figure

    co=Co(t_3Dplot,:);
    b1_co=xx_nucl(t_3Dplot,startB1Co:endB1Co);
    b2_co=xx_nucl(t_3Dplot,startB2Co:endB2Co);
    
    co_plot=zeros(n_bins);
    b1_co_plot=co_plot;
    
    for i1=1:n_bins
        for i2=1:n_bins
            index=ind(i2+n_bins*(i1-1));
            index_N=p.ind_scB1Co(i2+n_bins*(i1-1));
            if index >0 
                co_plot(i1,i2)=co(index);
                b1_co_plot(i1,i2)=b1_co(index);
                b2_co_plot(i1,i2)=b2_co(index);
            end
        end
    end

    b1_co_plot(co_plot<.1)=0;

    h=bar3((b1_co_plot)./co_plot);

    for k=1:length(h)
        h(k).LineStyle='none';
    end

    view([-150 30 90])

    axis_tick=p.age_bins;
    step_size_y=101;
    step_size_x=201;
    set(gca,'fontsize',28,'LineWidth',2,'ytick',...
        0:step_size_y-1:length(co),'yticklabel',...
        axis_tick(1):axis_tick(step_size_y):axis_tick(end)+Dtau,...
        'xtick',0:step_size_x-1:800,'xticklabel',...
        axis_tick(1):axis_tick(step_size_x):axis_tick(end)+Dtau,...
        'xlim',[0 601],'ylim',[0 601])%,'zlim',[0 1.2e6])

    if t_3Dplot ==t_plot(2)
        axis_tick=p.age_bins;
        step_size=51;
        set(gca,'fontsize',28,'LineWidth',2,'xtick',...
            0:step_size-1:length(co),'xticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'ytick',0:step_size-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'xlim',[0 251],'ylim',[0 251])
    else
        axis_tick=p.age_bins;
        step_size_y=101;
        step_size_x=201;
        set(gca,'fontsize',28,'LineWidth',2,'ytick',...
            0:step_size_y-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size_y):axis_tick(end)+Dtau,...
            'xtick',0:step_size_x-1:800,'xticklabel',...
            axis_tick(1):axis_tick(step_size_x):axis_tick(end)+Dtau,...
            'xlim',[0 601],'ylim',[0 601])
    end

end
set(gcf,'Position',[348   259   798   586])
colormap default

% Figure 7c
startB1Co=p.startB1Co;
endB1Co=p.endB1Co;
startB2Co=p.startB2Co;
endB2Co=p.endB2Co;
ind=p.ind;

for t_3Dplot=t_plot(4) 
    figure
    co=Co(t_3Dplot,:);
    b1_co=xx_nucl(t_3Dplot,startB1Co:endB1Co);
    b2_co=xx_nucl(t_3Dplot,startB2Co:endB2Co);
    co_plot=zeros(n_bins);
    b1_co_plot=co_plot;
    b2_co_plot=co_plot;
    
    for i1=1:n_bins
        for i2=1:n_bins
            index=ind(i2+n_bins*(i1-1));
            if index >0 
                co_plot(i1,i2)=co(index);
                b1_co_plot(i1,i2)=b1_co(index);
                b2_co_plot(i1,i2)=b2_co(index);
            end
        end
    end
    h=bar3(co_plot);
    view([-150 30 90])
    for k=1:length(h)
        h(k).LineStyle='none';
    end
    set(gca,'GridAlpha',0.1)
    if t_3Dplot ==t_plot(2)
        axis_tick=p.age_bins;
        step_size=51;
        set(gca,'fontsize',28,'LineWidth',2,'xtick',...
            0:step_size-1:length(co),'xticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'ytick',0:step_size-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'xlim',[0 251],'ylim',[0 251])
    else
        axis_tick=p.age_bins;
        step_size_y=101;
        step_size_x=201;
        set(gca,'fontsize',28,'LineWidth',2,'ytick',...
            0:step_size_y-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size_y):axis_tick(end)+Dtau,...
            'xtick',0:step_size_x-1:800,'xticklabel',...
            axis_tick(1):axis_tick(step_size_x):axis_tick(end)+Dtau,...
            'xlim',[0 601],'ylim',[0 601])
    end
    set(gcf,'Position',[348   259   798   586])
    colormap winter
end

% Figure 7d
startB1Co=p.startB1Co;
endB1Co=p.endB1Co;
startB2Co=p.startB2Co;
endB2Co=p.endB2Co;
ind=p.ind;

for t_3Dplot=t_plot(4)
    figure

    co=Co(t_3Dplot,:);
    b1_co=xx_nucl(t_3Dplot,startB1Co:endB1Co);
    b2_co=xx_nucl(t_3Dplot,startB2Co:endB2Co);
    
    co_plot=zeros(n_bins);
    b1_co_plot=co_plot;
    
    for i1=1:n_bins
        for i2=1:n_bins
            index=ind(i2+n_bins*(i1-1));
            index_N=p.ind_scB1Co(i2+n_bins*(i1-1));
            if index >0 
                co_plot(i1,i2)=co(index);
                b1_co_plot(i1,i2)=b1_co(index);
                b2_co_plot(i1,i2)=b2_co(index);
            end
        end
    end

    b1_co_plot(co_plot<.1)=0;

    h=bar3((b1_co_plot)./co_plot);

    for k=1:length(h)
        h(k).LineStyle='none';
    end

    view([-150 30 90])

    axis_tick=p.age_bins;
    step_size_y=101;
    step_size_x=201;
    set(gca,'fontsize',28,'LineWidth',2,'ytick',...
        0:step_size_y-1:length(co),'yticklabel',...
        axis_tick(1):axis_tick(step_size_y):axis_tick(end)+Dtau,...
        'xtick',0:step_size_x-1:800,'xticklabel',...
        axis_tick(1):axis_tick(step_size_x):axis_tick(end)+Dtau,...
        'xlim',[0 601],'ylim',[0 601])%,'zlim',[0 1.2e6])

    if t_3Dplot ==t_plot(2)
        axis_tick=p.age_bins;
        step_size=51;
        set(gca,'fontsize',28,'LineWidth',2,'xtick',...
            0:step_size-1:length(co),'xticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'ytick',0:step_size-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size):axis_tick(end),...
            'xlim',[0 251],'ylim',[0 251])
    else
        axis_tick=p.age_bins;
        step_size_y=101;
        step_size_x=201;
        set(gca,'fontsize',28,'LineWidth',2,'ytick',...
            0:step_size_y-1:length(co),'yticklabel',...
            axis_tick(1):axis_tick(step_size_y):axis_tick(end)+Dtau,...
            'xtick',0:step_size_x-1:800,'xticklabel',...
            axis_tick(1):axis_tick(step_size_x):axis_tick(end)+Dtau,...
            'xlim',[0 601],'ylim',[0 601])
    end

end
set(gcf,'Position',[348   259   798   586])
colormap winter