function p = initialization_parameters(age_viab,age_end_inf,Dtau,Drt,C0,MOI_1,MOI_2,scaling)
    
    % binding decay fitting parameters
    k_bind0_I1=6.3e-7*scaling; % mL/cell/h 
    tau0_bind_I1=1.80;
    beta_bind_I1=0.5;
    
    % cell growth and death 
    mu_T=0.028; % 1/h
    mu_I1=0;
    k_deathT=8e-5; % 1/h
    tau_death_I1=24;
    k_death_I1=0.0029; % 1/h
    
    % virus degradation in medium
    k_degrBV_1=7e-3; %1/h
    
    % BV internalization
    k_int_I1=0.01*60; %1/h
    eta_I1=0.5;
    
    % replication
    k_repl_I1=.7318;  % 1/h
    replStart_I1=6; % h
    replEnd_I1=18; % h
    
    % budding
    k_rel_I1=5;   
    tau_rel_on_I1=18;
    tau_rel_off_I1=72;

    % viral genome degradation in nucleus
    k_d_N1=0; %1/h

    % Substrate
    K_s=1.3*1e3/scaling; % nmol/mL=1.3 mM raghunand and dale, 1996
    Ys_T=1.2e-4; %nmol/cell/h raghunand and dale, 1999
    Ys_I1=0; 
    
    %% Mesh and state vectors indexing
    % Mesh
    n_bins=age_viab/Dtau;
    n_bins_inf=age_end_inf/Dtau;
    age_bins=0:Dtau:age_viab-Dtau;
    
    % first state vector
    startI1=3;
    endI1=2+n_bins;
    
    % other state vectors: internalized virus
    startB1=1;
    endB1=n_bins;
    
    ind=zeros(1,n_bins);
    
    %% Parameters infected cells
    f_rel_I1=zeros(1,length(age_bins));
    f_rel_I1(age_bins>=tau_rel_on_I1)=1;
    f_rel_I1(age_bins>=tau_rel_off_I1)=0;
    
    f_repl_I1=(0.5+0.5*tanh((age_bins-replStart_I1)/0.3)).*(0.5+0.5*tanh((replEnd_I1-age_bins)/1)); % smooth transition
    
    n1=0.5+0.5*(tanh((age_bins(1:n_bins_inf)-tau0_bind_I1)/.01));
    k_bind_I1=k_bind0_I1*((1-n1)+n1.*exp(-beta_bind_I1*(age_bins(1:n_bins_inf)-tau0_bind_I1)));
    
    % correction to reallocation when age_end_inf = age_viability
    if age_end_inf == age_viab
        a=-1;
    else
        a=0;
    end

    %% Initialization of solver and state vector
    hmax0 = Dtau; % CFL condition % ode23 default: 0.1*(T-t0);
    rtol = 1.e-10;
    atol = 1.e-6;
    threshold=atol/rtol;
    
    %% Prepare parameters object to pass to solver
    p.startI1=startI1;
    p.endI1=endI1;

    p.startB1=startB1;

    p.endB1=endB1;
    p.age_end_inf=age_end_inf;
    p.a=a;
    
    p.ind=ind;
    p.n_bins=n_bins;
    
    p.n_bins_inf=n_bins_inf;
    p.age_bins=age_bins;
    
    p.tau_death_I1=tau_death_I1;
    p.k_death_I1=k_death_I1;
    p.k_bind_I1=k_bind_I1;
    
    p.mu_T=mu_T;
    p.mu_I1=mu_I1;
    p.k_deathT=k_deathT;
    p.k_int_I1=k_int_I1;
    p.eta_I1=eta_I1;
    p.f_repl_I1=f_repl_I1;
    p.k_repl_I1=k_repl_I1;
    p.f_rel_I1=f_rel_I1;
    p.k_rel_I1=k_rel_I1;
    p.k_degrBV_1=k_degrBV_1;

    p.tau_rel_on_I1=tau_rel_on_I1;
    p.tau_rel_off_I1=tau_rel_off_I1;
    p.replStart_I1=replStart_I1;
    p.replEnd_I1=replEnd_I1;
    p.tau0_bind_I1=tau0_bind_I1;
    p.k_bind0_I1=k_bind0_I1;
    p.beta_bind_I1=beta_bind_I1;

    p.k_d_N1=k_d_N1;
    
    p.K_s=K_s;
    p.Ys_T=Ys_T;
    p.Ys_I1=Ys_I1; 
    
    p.threshold=threshold;
    p.rtol=rtol;
    p.Dtau=Dtau;
    p.hmax0=hmax0;