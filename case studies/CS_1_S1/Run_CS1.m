% Case study 1
clear, close all

MOIs=[1 0.1 0.001];

ylims_I=[5e4 5e4 5000 5000];
cols=[{[189 148 196]/255}  {[51 74 159]/255} {[163 202 235]/255}  {[0 169 192]/255}];
lines=[{'-'},{'-.'},{':'}];
lws=[3 4 5];

% numeric_scheme='RK45';
numeric_scheme='RK23';

%% Generate Figure 6
for kk=1:length(MOIs)
    lw=lws(kk);
    col=cols{kk};
    line=lines{kk};
    MOI_1=MOIs(kk);
    CS1_plots(lw,col,line,MOI_1,ylims_I,numeric_scheme);
end
% save_plots();

function CS1_plots(lw,col,line,MOI_1,ylims_I,numeric_scheme)

    scaling=1e7;    
    pos=[2079   519   473   245];
    
    % Simulation parameters
    t_final=24*10;
    C0=1.5e6;
    C0_nv=0;
    S0=15*1e3; % nmol/mL
    Cin=3e6;
    Din=0.01;
    Sin=15*1e3/scaling; % nmol/mL
    MOI_2=0;
    age_end_inf=20;
    age_viab=140;    
    Dt = 1; % h: control interval
    
    %% Control law
    Dtau=.1; % mesh grid
    tt=0:Dt:t_final;
    
    DD=ones(1,length(tt))*Din;
    CCin=ones(1,length(tt))*Cin;
    r_bleed_vect=ones(1,length(tt));
        
    %% Initialization
    p=def_parameters(age_viab,age_end_inf,Dtau,Dt,C0,MOI_1,MOI_2,scaling);
    
    n_bins=p.n_bins;
    t=0;
    V0_1=C0*MOI_1; % #/mL
    x0=[C0 V0_1 zeros(1,n_bins) C0_nv S0]/scaling;
    x0_bind=zeros(1,n_bins);
    x0_nucl=x0_bind;
    sum_h=0;
    
    % Preallocate solution vectors
    tt=0:Dt:t_final;
    xx=zeros(length(tt),2+n_bins+2);
    xx_bind=zeros(length(tt),p.endB1);
    xx_nucl=xx_bind;
    xx(1,:)=x0*scaling;
        
    %% Simulation
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
                x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling);
    
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
    
    %% Plots
    I1=xx(:,p.startI1:p.endI1);

    figure(1)
    hold on
    plot(tt/24,sum(I1,2),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
    set(gcf,'Position',pos)
    xlabel('Time [d]')
    ylabel('Infected [cell/mL]')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(2)
    hold on
    plot(tt/24,xx(:,2),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22,'ylim',[0 2e8]) %,'xticklabel',[])
    set(gcf,'Position',pos)
    ylabel('Virion [PFU/mL]')
    xlabel('Time [d]')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1],'xtick',0:2:20)
    box on
    
    figure(3)
    hold on
    plot(tt/24,xx(:,end)/1e3,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
    set(gcf,'Position',pos)
    ylabel('Glucose [mM]')
    xlabel('Time [d]')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(4)
    hold on 
    plot(tt/24,xx(:,1),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'ylim',[0 4e6]) %,'xticklabel',[])
    set(gcf,'Position',pos)
    ylabel('Uninfected [cell/mL]')
    xlabel('Time [d]')
    box on
    
    DNA1=xx_nucl(:,p.startB1:p.endB1);

    t_plot=[48 96 144 240]+1;
    pos2=[ 323   454   560   150];

    % infection age distribution
    for i=1:length(t_plot)
        figure(5+i)
        hold on
        ylim_I=ylims_I(i);
        plot(p.age_bins,I1(t_plot(i),:),'linewidth',lw,'Color',col,'LineStyle',line)
        set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],'ylim',[0 ylim_I],'fontsize',22,'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
        box on
        set(gcf,'Position',pos2);
    end
    set(gca,'linewidth',2,'xticklabel',0:20:140,'fontsize',22,'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
        
    % DNA/cell
    for i=1:length(t_plot) 
        figure(9+i)
        hold on
        vDNA_plot=DNA1(t_plot(i),:)./(I1(t_plot(i),:));
        age_bins_plot=p.age_bins;
        age_bins_plot(vDNA_plot<.1)=[];
        vDNA_plot(vDNA_plot<.1)=[];
        plot(age_bins_plot,vDNA_plot,'linewidth',lw,'Color',col,'LineStyle',line)
        set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],'fontsize',22,'ylim',[0 4e5],'xlim',[0 100],'YScale','log') %,'xticklabel',[])
        set(gcf,'Position',pos2);
        box on
        set(gca,'ylim',[1 1e9],'ytick',[1 1e5])
    end
    set(gca,'xticklabel',0:20:140) 
end

function save_plots()
    figure(1)
    saveas(gcf,'I_CS1.png')
    
    figure(2)
    saveas(gcf,'V_CS1.png')
    
    figure(3)
    set(gcf,'Position',[2922    520   473  229])
    saveas(gcf,'S_CS1.png')
    
    figure(4)
    saveas(gcf,'T_CS1.png')
    
    figure(6)
    set(gca,'ytick',[0 4e4])
    text(55, 50000*.8,'2 d post inoculation','FontSize',20)
    saveas(gcf,'age_CS1_1.png')
    
    figure(7)
    set(gca,'ytick',[0 4e4])
    text(55, 50000*.8,'4 d post inoculation','FontSize',20)
    saveas(gcf,'age_CS1_2.png')
    
    figure(8)
    set(gcf,'Position',[ 729   458   565   122])
    set(gca,'ytick',[0 4000])
    text(55, 5000*.8,'6 d post inoculation','FontSize',20)
    saveas(gcf,'age_CS1_3.png')
    
    figure(9)
    set(gcf,'Position',[  861   202   565   142])
    set(gca,'ytick',[0 4000])
    text(55, 5000*.8,'10 d post inoculation','FontSize',20)
    saveas(gcf,'age_CS1_4.png')
    
    figure(10)
    text(55, 1e9*.02,'2 d post inoculation','FontSize',20)
    set(gcf,'Position',[ 1608         123         560         120])
    saveas(gcf,'DNA_pc_CS1_1.png')
    
    figure(11)
    text(55, 1e9*.02,'4 d post inoculation','FontSize',20)
    set(gcf,'Position',[ 1608         123         560         120])
    saveas(gcf,'DNA_pc_CS1_2.png')
    
    figure(12)
    text(55,1e9*.02,'6 d post inoculation','FontSize',20)
    set(gcf,'Position',[ 1608         123         560         120])
    saveas(gcf,'DNA_pc_CS1_3.png')
    
    figure(13)
    text(55, 1e9*.02,'10 d post inoculation','FontSize',20)
    set(gcf,'Position',[ 1703        -308         560         139])
    saveas(gcf,'DNA_pc_CS1_4.png')
end