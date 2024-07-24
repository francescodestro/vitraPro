%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function for defining the model parameters 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Detailed description of the parameters: Table S2 in Destro and Braatz, 2024 
%       For STV/DIP systems: STV = virus 1, DIP = virus 2.
%       For systems with two STVs, no predefined order.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = def_parameters(age_viab,age_end_inf,Dtau,scaling)   

    % cell growth and death 
    mu_T=0.028; % 1/h
    mu_I1=0; % 1/h
    k_deathT=8e-5; % 1/h
    tau_death_I1=24; % hpi
    k_death_I1=0.0029; % 1/h

    % viral binding 
    k_bind0_I1=8e-7*scaling; % mL/cell/h 
    tau0_bind_I1=1.80; % hpi
    beta_bind_I1=.5; % –
    
    % virus internalization
    k_int_I1=0.01*60; %1/h
    eta_I1=.5; % -
    
    % viral replication
    k_repl_I1=.7318;  % 1/h
    replStart_I1=6; % hpi
    replEnd_I1=18; % hpi

    % viral progeny release
    k_rel_I1=5;   % PFU/cell/h
    tau_rel_on_I1=18; % hpi
    tau_rel_off_I1=72; % hpi

    % viral genome degradation in nucleus
    k_d_N1=0; % 1/h

    % viral degradation in supernatant
    k_d_V1=7e-3; %1/h

    % Substrate
    K_s=1.3*1e3; % nmol/mL
    Ys_T=1.2e-4; % nmol/cell/h 
    Ys_I1=0;  % nmol/cell/h 

    % Nonviable cell lysis
    k_lys=0; % 1/h

    %% Don't edit part below
    % Mesh and state vectors indexing
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
    
    K_s=K_s/scaling;
    
    %% Parameters infected cells
    f_rel_I1=zeros(1,length(age_bins));
    f_rel_I1(age_bins>=tau_rel_on_I1)=1;
    f_rel_I1(age_bins>=tau_rel_off_I1)=0;
    f_rel_I1(end)=0;

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
    hmax0 = Dtau; % CFL condition 
    rtol = 1.e-6;
    atol = 1.e-6;
    threshold=atol/rtol;
    
    %% Prepare parameters object to pass to solver
    p.startI1=startI1;
    p.endI1=endI1;
    p.startB1=startB1;
    p.endB1=endB1;
    p.age_end_inf=age_end_inf;
    p.a=a;
    
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
  
    p.k_degrBV_1=k_d_V1;

    p.k_lys=k_lys;
    p.k_d_N1=k_d_N1;

    p.K_s=K_s;
    p.Ys_T=Ys_T;
    p.Ys_I1=Ys_I1; 

    p.threshold=threshold;
    p.rtol=rtol;
    p.Dtau=Dtau;
    p.hmax0=hmax0;

end
