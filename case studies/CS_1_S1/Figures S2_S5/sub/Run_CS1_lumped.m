function Run_CS2_lumped(par)

    % Case study 2
    if nargin<1
        par=exp([-1.5763   -2.7799    1.7034]);
    end

    lw=4.5;
    % col=[0.8500, 0.3250, 0.0980];
    col=[0.4660, 0.6740, 0.1880]*1.2;
    line=':';
    
    scaling=1e8;
    
    t_final=24*10;
    C0=1.5e6;
    C0_nv=0;
    S0=15*1e3; % nmol/mL
    Cin=3e6;
    Din=0.01;
    Sin=15*1e3/scaling; % nmol/mL
    MOI_1=1;
    MOI_2=0;
    age_end_inf=20;
    age_viab=140;    
    Dt = .1; % h: control interval
    
    %% Control law
    Dtau=.1; % mesh grid
    tt=0:Dt:t_final;
    
    DD=ones(1,length(tt))*Din;
    CCin=ones(1,length(tt))*Cin;
    r_bleed_vect=ones(1,length(tt));
    
    %% Initialization
    p=initialization_parameters(age_viab,age_end_inf,Dtau,Dt,C0,MOI_1,MOI_2,scaling);
    
    t=0;
    V0_1=C0*MOI_1; % #/mL
    x0=[C0 V0_1 0 S0]/scaling;
    
    % Preallocate solution vectors
    tt=0:Dt:t_final;
    xx=zeros(length(tt),4);
    xx(1,:)=x0*scaling;
    
    %% Parameters lumped model
    
    % par=exp([-1.5763   -2.7799    1.7034]);
    
    mu=p.mu_T;
    k_bind=p.k_bind0_I1;
    k_deathT=p.k_deathT;
    k_deathI=p.k_death_I1;
    k_rel=p.k_rel_I1;
    K_s=p.K_s;
    Ys_T=p.Ys_T;
    Ys_I=p.Ys_I1;
    k_degrV=p.k_degrBV_1;
    
    k_bind=par(1);
    k_deathI=par(2);
    k_rel=par(3);
    
    tic
    %% Simulation
    for i = 2:length(tt)
    
        D=DD(i);
        Cin=CCin(i);
        r_bleed=r_bleed_vect(i);
    
        [t_out,x] = ode23s(@inf_model_lumped,[t t+Dt],x0,[],mu,k_bind,k_deathT,...
            k_deathI,k_rel,K_s, Ys_T, Ys_I, k_degrV, Cin, r_bleed, Sin, D, scaling);
    
        % initialize following step
        x0=x(end,:);
        
        t=t_out(end);
    
        % save new results
        index=i;
        xx(index,:)=x(end,:)*scaling;
    
    end
    toc
    
    %% Plots
    figure(1)
    hold on
    plot(tt/24,xx(:,3),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
    xlabel('Time [d]')
    ylabel('Infected [cell/mL]')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(2)
    hold on
    plot(tt/24,xx(:,2),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22,'ylim',[0 8e8]) %,'xticklabel',[])
    ylabel('Virion [PFU/mL]')
    xlabel('Time [d]')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(3)
    hold on
    plot(tt/24,xx(:,end)/1e3,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
    ylabel('Glucose [mM]')
    xlabel('Time [d]')
    set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
    box on
    
    figure(4)
    hold on 
    plot(tt/24,xx(:,1),'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'linewidth',2,'fontsize',22) %,'ylim',[0 4e6]) %,'xticklabel',[])
    ylabel('Uninfected [cell/mL]')
    xlabel('Time [d]')
    box on  
