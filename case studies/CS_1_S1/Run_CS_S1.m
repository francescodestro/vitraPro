% Case study S1
clear

MOIs=[1 0.1 0.001];

cols=[{[189 148 196]/255} {[51 74 159]/255} {[163 202 235]/255} {[0 169 192]/255} {[0 169 192]/255}];
lines=[{'-'},{'-.'},{':'}];
lws=[3 4 5]; 

% numeric_scheme='RK45';
numeric_scheme='RK23';

%% Generate Figure S1
for kk=1:length(MOIs)
    CS_S1_plots(kk,MOIs,cols,lines,lws,numeric_scheme);
end
% save_plots();

%% % Generate Figure 4
% for kk=1
%     parameters_plot(kk,MOIs,cols,lines,lws);
% end

function CS_S1_plots(kk,MOIs,cols,lines,lws,numeric_scheme)
    col=cols{kk};
    line=lines{kk};
    MOI_1=MOIs(kk);
    lw=lws(kk);
    scaling=1e6;

    pos=[2079   519     473     245];
    t_final=24*10;
    C0=1.5e6;
    C0_nv=0;
    S0=15*1e3; % nmol/mL
    Cin=3e6;
    Din=0.0;
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
        
    tic
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
    %     toc
    end
    toc
    
    %% Plots
    I1=xx(:,p.startI1:p.endI1);

    figure(1)
    hold on
    plot(tt/24,sum(I1,2),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
    set(gcf,'Position',pos)
    xlabel('Time [d]','FontName','arial')
    ylabel('Infected [cell/mL]','FontName','arial')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(2)
    hold on
    plot(tt/24,xx(:,2),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22,'ylim',[0 8e8]) %,'xticklabel',[])
    set(gcf,'Position',pos)
    ylabel('Virion [PFU/mL]','FontName','arial')
    xlabel('Time [d]','FontName','arial')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(3)
    hold on
    plot(tt/24,xx(:,end)/1e3,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
    set(gcf,'Position',pos)
    ylabel('Glucose [mM]','FontName','arial')
    xlabel('Time [d]','FontName','arial')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(4)
    hold on 
    plot(tt/24,xx(:,1),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'ylim',[0 4e6]) %,'xticklabel',[])
    set(gcf,'Position',pos)
    ylabel('Uninfected [cell/mL]','FontName','arial')
    xlabel('Time [d]','FontName','arial')
    box on
    
    DNA1=xx_nucl(:,p.startB1:p.endB1);

    t_plot=[25 49 73 97];
    pos2=[ 323   454   560   150];

    % infection age distribution
    for i=1:length(t_plot)
        figure(5+i)
        hold on
        I1_plot=I1(t_plot(i),:);
        plot(p.age_bins,I1_plot,'linewidth',lw,'Color',col,'LineStyle',line)
        set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],...
            'fontsize',22,'ylim',[0 1.5e5],'xlim',[0 100],'ytick',[0 5 10]*1e4,...
            'YScale','linear','FontName','arial') %,'xticklabel',[])
        set(gcf,'Position',pos2);
        box on
    end
    set(gca,'linewidth',2,'xticklabel',0:20:140,'fontsize',22,'ylim',[0 1.5e5],'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
    
    % DNA/cell
    for i=1:length(t_plot) 
        figure(9+i)
        hold on
        vDNA_plot=DNA1(t_plot(i),:)./(I1(t_plot(i),:));
        age_bins_plot=p.age_bins;
%         age_bins_plot(vDNA_plot<.1)=[];
%         vDNA_plot(vDNA_plot<.1)=[];
        age_bins_plot(I1(t_plot(i),:)<10)=0;
        vDNA_plot(I1(t_plot(i),:)<10)=0;

        plot(age_bins_plot(),vDNA_plot(),'linewidth',lw,'Color',col,'LineStyle',line)
        set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],...
            'fontsize',22,'xlim',[0 100],'YScale','log') %,'xticklabel',[])
        set(gcf,'Position',pos2);
        box on 
        set(gca,'ylim',[1 1e9],'ytick',[1 1e5],'FontName','arial')
%         set(gca,'ylim',[0 4e5],'YScale','linear')
    end
    set(gca,'xticklabel',0:20:140) 
end

function save_plots()
    figure(1)
    set(gca,'FontName','arial')
    saveas(gcf,'I_CS_S1.png')
    
    figure(2)
    set(gca,'FontName','arial')
    saveas(gcf,'V_CS_S1.png')
    
    figure(3)
    set(gcf,'Position',[2922    520   473  229])
    set(gca,'FontName','arial')
    saveas(gcf,'S_CS_S1.png')
    
    figure(4)
    set(gca,'FontName','arial')
    saveas(gcf,'T_CS_S1.png')
    
    figure(6)
    set(gca,'ytick',[0 5e4 10e4])
    set(gca,'ylim',[0 15e4])
    text(55, 15e4*.8,'1 d post inoculation','FontSize',20,'FontName','arial')
    set(gca,'FontName','arial')
    saveas(gcf,'age_CS_S1_1.png')
    
    figure(7)
    set(gca,'ytick',[0 5e4])
    set(gca,'ylim',[0 7.5e4])
    text(55, 7.5e4*.8,'2 d post inoculation','FontSize',20,'FontName','arial')
    set(gca,'FontName','arial')
    saveas(gcf,'age_CS_S1_2.png')
    
    figure(8)
    set(gca,'ylim',[0 7.5e4])
    set(gca,'ytick',[0 5e4])
    text(55, 7.5e4*.8,'3 d post inoculation','FontSize',20,'FontName','arial')
    set(gca,'FontName','arial')
    saveas(gcf,'age_CS_S1_3.png')
    
    figure(9)
    set(gcf,'Position',[743   368   560   162])
    set(gca,'ylim',[0 7.5e4])
    set(gca,'ytick',[0 5e4])
    text(55, 7.5e4*.8,'4 d post inoculation','FontSize',20,'FontName','arial')
    set(gca,'FontName','arial')
    saveas(gcf,'age_CS_S1_4.png')
    
    figure(10)
    text(55, 1e9*.02,'1 d post inoculation','FontSize',20,'FontName','arial')
    set(gcf,'Position',[ 1608         123         560         123])
    set(gca,'FontName','arial')
    saveas(gcf,'DNA_pc_CS_S1_1.png')
    
    figure(11)
    text(55, 1e9*.02,'2 d post inoculation','FontSize',20,'FontName','arial')
    set(gcf,'Position',[ 1608         123         560         123])
    set(gca,'FontName','arial')
    saveas(gcf,'DNA_pc_CS_S1_2.png')
    
    figure(12)
    text(55,1e9*.02,'3 d post inoculation','FontSize',20,'FontName','arial')
    set(gcf,'Position',[ 1608         123         560         123])
    set(gca,'FontName','arial')
    saveas(gcf,'DNA_pc_CS_S1_3.png')
    
    figure(13)
    text(55, 1e9*.02,'4 d post inoculation','FontSize',20,'FontName','arial')
    set(gca,'FontName','arial')
    set(gcf,'Position',[ 1703        -308         560         142])
%     set(gcf,'Position',[743   368   560   162])
    saveas(gcf,'DNA_pc_CS_S1_4.png')

end

function parameters_plot(kk,MOIs,cols,lines,lws)

    col=cols{5};
    line=lines{kk};
    MOI_1=MOIs(kk);
    lw=lws(kk);
    scaling=1e8;

    % Simulation parameters
    C0=1.5e6;

    MOI_2=0;
    age_end_inf=20;
    age_viab=140;    
    Dt = 1; % h: control interval
    
    Dtau=.1; % mesh grid
        
    p=initialization_parameters(age_viab,age_end_inf,Dtau,Dt,C0,MOI_1,MOI_2,scaling);
    
    pos2=[ 323   454   500   250];

    figure
    k_bindI1=zeros(1,length(p.age_bins));
    k_bindI1(1:p.n_bins_inf)=p.k_bind_I1;
    plot(p.age_bins,k_bindI1/scaling,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'YAxisLocation','left','linewidth',2,'fontsize',22,...
        'ylim',[0 8e-7],'xlim',[0 100],'XMinorTick','on','FontName','arial')%,'ylim',[0 2.5e9],'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
    set(gcf,'Position',pos2);
    xlabel('Infection age [hpi]','FontName','arial')
    ylabel('k̅_{b,I} [mL cell^{-1} h^{-1}]','Interpreter','tex')
    box on

    figure
    plot(p.age_bins,p.k_rel_I1.*p.f_rel_I1,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'YAxisLocation','left','linewidth',2,'fontsize',22,...
        'ylim',[0 8],'xlim',[0 100],'XMinorTick','on','FontName','arial')%,'ylim',[0 2.5e9],'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
    set(gcf,'Position',pos2);
    xlabel('Infection age [hpi]','FontName','arial')
    ylabel('k̅_{v} [PFU cell^{-1} h^{-1}]','Interpreter','tex')
    box on

    figure
    plot(p.age_bins,p.k_repl_I1.*p.f_repl_I1,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'YAxisLocation','left','linewidth',2,'fontsize',22,'ylim',[0 .9],...
        'ytick',0:.3:1,'xlim',[0 100],'XMinorTick','on','FontName','arial')%,'ylim',[0 2.5e9],'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
    set(gcf,'Position',pos2);
    xlabel('Infection age [hpi]','FontName','arial')
    ylabel('k̅_{r} [h^{-1}]','Interpreter','tex')
    box on

    figure
    DNA=linspace(1,5e5,500);
    k_deathI1=p.k_death_I1*log(DNA);
    semilogx(DNA,k_deathI1,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'YAxisLocation','left','linewidth',2,'fontsize',22,'ylim',[0 .06],...
        'ytick',0:.02:1,'xlim',[0 5e5],'XMinorTick','on','FontName','arial')%,'ylim',[0 2.5e9],'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
    set(gcf,'Position',pos2);
    xlabel('Viral genome [# cell^{-1}]','FontName','arial')
    ylabel('k̅_{d,I} [h^{-1}]','Interpreter','tex')
    box on
end