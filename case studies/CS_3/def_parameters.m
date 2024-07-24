function p = def_parameters(age_viab,age_end_inf,Dtau,scaling)   

    % binding decay fitting parameters
    k_bind0_I1=1.5e-6*scaling;
    tau0_bind_I1=1.5;
    beta_bind_I1=100;
    
    k_bind0_I2=k_bind0_I1;
    tau0_bind_I2=15;
    beta_bind_I2=10;
    
    tau0_bind_Co1=tau0_bind_I1;
    beta_bind_Co1=beta_bind_I1;
    
    tau0_bind_Co2=tau0_bind_I1;
    beta_bind_Co2=beta_bind_I1;
    
    % cell growth, death 
    mu_T= 0.02;
    mu_I1=0;
    mu_I2=mu_T;
    mu_Co=0;
    k_deathT= 0.0086;
    tau_death_I1=6.5;
    tau_death_I2=tau_death_I1;
    tau_death_Co=tau_death_I1;
    k_death_I1=0.0300;
    k_death_I2=k_deathT;
    k_death_Co=k_death_I1;
    f_death_Co1=1; 
    f_death_Co2=0;
    
    % virus degradation in medium
    k_degrBV_1=0.17;
    k_degrBV_2=k_degrBV_1; 
    
    % Viral internalization
    k_int_I1=3; 
    k_int_I2=k_int_I1;
    k_int_Co1=k_int_I1; 
    k_int_Co2=k_int_I1; 
    eta_I1=.5;
    eta_I2=eta_I1;
    eta_Co1=eta_I1;
    eta_Co2=eta_I1;
    
    % replication
    k_repl_I1=1.7839; 
    replStart_I1=2.2158;
    replEnd_I1=6.0794;
    k_repl_I2=0; 
    replStart_I2=6; 
    replEnd_I2=18; 
    k_repl_Co1=k_repl_I1;  
    replStart_Co1=replStart_I1; 
    replEnd_Co1=replEnd_I1; 
    k_repl_Co2=k_repl_I1*1.3132; 
    replStart_Co2=replStart_I1;
    replEnd_Co2=replEnd_I1; 
    
    % budding
    k_rel_I1=32.9108;   
    tau_rel_on_I1=6.5;
    tau_rel_off_I1=60;
    k_rel_I2=0;   
    tau_rel_on_I2=10;
    tau_rel_off_I2=72;
    k_rel_Co=k_rel_I1;   
    tau_rel_on_Co=tau_rel_on_I1;
    tau_rel_off_Co=tau_rel_off_I1;
    randV2_I1=0;
    randV1_I2=0;

    % viral genome degradation in nucleus
    k_d_N1=0; 
    k_d_N2=0; 

    % nonviable cell lysis
    k_lys=0.15;

    % Substrate
    K_s=0;
    Ys_T=0;
    Ys_I1=0; 
    Ys_I2=0; 
    Ys_C=0;

    %% Mesh and state vectors indexing
    % Mesh
    n_bins=age_viab/Dtau;
    n_bins_inf=age_end_inf/Dtau;
    age_bins=0:Dtau:age_viab-Dtau;
    
    % first state vector
    startI1=4;
    endI1=3+n_bins;
    startI2=4+n_bins;
    endI2=3+2*n_bins;
    startCo=4+2*n_bins;
    endCo=3+2*n_bins+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);
    
    % other state vectors: internalized virus
    startB1=1;
    endB1=n_bins;
    startB2=1+n_bins;
    endB2=2*n_bins;
    startB1Co=2*n_bins+1;
    endB1Co=2*n_bins+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);
    startB2Co=endB1Co+1;
    endB2Co=endB1Co+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);
    
    ind=zeros(1,n_bins^2);
    ind_scCo=ind;
    ind_scB1Co=ind;
    ind_scB2Co=ind;
    
    %% Parameters for coinfected cells
    k_bind0_Co1=k_bind0_I1;
    k_bind0_Co2=k_bind0_I1;

    k_bind_Co1=ind;
    k_bind_Co2=ind;
    
    age_co=ind;
    f_repl_Co1=ind;
    f_repl_Co2=ind;
    
    % indexes co-infected cells bins and binding constant of (co-)infected cells
    i=1;
    for i1=1:n_bins
        for i2=1:n_bins
            if (i1-i2)^2 < n_bins_inf^2
                ind(i2+n_bins*(i1-1))=i; % entry of dCodt, age_co, f_repl, f_rel
                ind_scCo(i2+n_bins*(i1-1))=i+startCo-1; % entry of dxdt
                ind_scB1Co(i2+n_bins*(i1-1))=i+startB1Co-1; % entry of dbind_dt and dnucl_dt
                ind_scB2Co(i2+n_bins*(i1-1))=i+startB2Co-1; % entry of dbind_dt and dnucl_dt
    
                current_age=((i1))*Dtau;
    
                if i2==n_bins
                    current_age=((i2)-1)*Dtau;
                end
    
                n1=0.5+0.5*(tanh((current_age-tau0_bind_Co1)/.01));
                k_bind_Co1(i2+n_bins*(i1-1))=k_bind0_Co1*((1-n1)+n1.*exp(-beta_bind_Co1*(current_age-tau0_bind_Co1)));
    
                n2=0.5+0.5*(tanh((current_age-tau0_bind_Co2)/.01));
                k_bind_Co2(i2+n_bins*(i1-1))=k_bind0_Co2*((1-n1)+n2.*exp(-beta_bind_Co2*(current_age-tau0_bind_Co2)));
    
                f_repl_Co1(i2+n_bins*(i1-1))=(0.5+0.5*tanh((current_age-replStart_Co1)/0.3)).*...
                    (0.5+0.5*tanh((replEnd_Co1-current_age)/1));
                f_repl_Co2(i2+n_bins*(i1-1))=(0.5+0.5*tanh((current_age-replStart_Co2)/0.3)).*...
                    (0.5+0.5*tanh((replEnd_Co2-current_age)/1));
    
                age_co(i2+n_bins*(i1-1))=current_age;
                i=i+1;
            end
        end
    end
    
    age_co=age_co(ind>0);
    f_repl_Co1=f_repl_Co1(ind>0);
    f_repl_Co2=f_repl_Co2(ind>0);
    k_bind_Co1=k_bind_Co1(ind>0);
    k_bind_Co2=k_bind_Co2(ind>0);
    
    f_rel_Co=zeros(1,length(age_co));
    f_rel_Co(age_co>=tau_rel_on_Co)=1;
    f_rel_Co(age_co<=tau_rel_on_Co)=0;
    f_rel_Co(age_co>=tau_rel_off_Co)=0;
    f_rel_Co(age_co==(age_viab-Dtau))=0;
    
    %% Parameters infected cells
    f_rel_I1=zeros(1,length(age_bins));
    f_rel_I1(age_bins>=tau_rel_on_I1)=1;
    f_rel_I1(age_bins>=tau_rel_off_I1)=0;
    f_rel_I1(end)=0;

    f_rel_I2=zeros(1,length(age_bins));
    f_rel_I2(age_bins>=tau_rel_on_I2)=1;
    f_rel_I2(age_bins>=tau_rel_off_I2)=0;
    f_rel_I2(end)=0;
    
    f_repl_I1=(0.5+0.5*tanh((age_bins-replStart_I1)/0.3)).*(0.5+0.5*tanh((replEnd_I1-age_bins)/1)); % smooth transition
    f_repl_I2=(0.5+0.5*tanh((age_bins-replStart_I2)/0.3)).*(0.5+0.5*tanh((replEnd_I2-age_bins)/1)); % smooth transition
    
    n1=0.5+0.5*(tanh((age_bins(1:n_bins_inf)-tau0_bind_I1)/.01));
    k_bind_I1=k_bind0_I1*((1-n1)+n1.*exp(-beta_bind_I1*(age_bins(1:n_bins_inf)-tau0_bind_I1)));
    
    n2=0.5+0.5*(tanh((age_bins(1:n_bins_inf)-tau0_bind_I2)/.01));
    k_bind_I2=k_bind0_I2*((1-n2)+n2.*exp(-beta_bind_I2*(age_bins(1:n_bins_inf)-tau0_bind_I2)));
    
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
    p.startI2=startI2;
    p.endI2=endI2;
    p.startCo=startCo;
    p.endCo=endCo;
    p.startB1=startB1;
    p.startB2=startB2;
    p.endB1=endB1;
    p.endB2=endB2;
    p.endB1Co=endB1Co;
    p.endB2Co=endB2Co;
    p.startB1Co=startB1Co;
    p.startB2Co=startB2Co;
    p.ind_scCo=ind_scCo;
    p.age_end_inf=age_end_inf;
    p.a=a;
    
    p.ind=ind;
    p.ind_scB1Co=ind_scB1Co;
    p.ind_scB2Co=ind_scB2Co;
    
    p.n_bins=n_bins;
    
    p.n_bins_inf=n_bins_inf;
    p.age_bins=age_bins;
    p.age_co=age_co;
    
    p.tau_death_I1=tau_death_I1;
    p.k_death_I1=k_death_I1;
    p.tau_death_I2=tau_death_I2;
    p.k_death_I2=k_death_I2;
    p.tau_death_Co=tau_death_Co;
    p.k_death_Co=k_death_Co;
    p.f_death_Co1=f_death_Co1;
    p.f_death_Co2=f_death_Co2;
    
    p.k_bind_I1=k_bind_I1;
    p.k_bind_I2=k_bind_I2;
    p.k_bind_Co1=k_bind_Co1;
    p.k_bind_Co2=k_bind_Co2;
    
    p.mu_T=mu_T;
    p.mu_I1=mu_I1;
    p.mu_I2=mu_I2;
    p.mu_Co=mu_Co;
    p.k_deathT=k_deathT;
    
    p.k_int_I1=k_int_I1;
    p.k_int_I2=k_int_I2;
    p.k_int_Co1=k_int_Co1;
    p.k_int_Co2=k_int_Co2;
    
    p.eta_I1=eta_I1;
    p.eta_I2=eta_I2;
    p.eta_Co1=eta_Co1;
    p.eta_Co2=eta_Co2;
    
    p.f_repl_I1=f_repl_I1;
    p.f_repl_I2=f_repl_I2;
    p.f_repl_Co1=f_repl_Co1;
    p.f_repl_Co2=f_repl_Co2;
    p.k_repl_I1=k_repl_I1;
    p.k_repl_I2=k_repl_I2;
    p.k_repl_Co1=k_repl_Co1;
    p.k_repl_Co2=k_repl_Co2;
    
    p.f_rel_I1=f_rel_I1;
    p.f_rel_I2=f_rel_I2;
    p.f_rel_Co=f_rel_Co;
    p.k_rel_I1=k_rel_I1;
    p.k_rel_I2=k_rel_I2;
    p.k_rel_Co=k_rel_Co;
    
    p.k_degrBV_1=k_degrBV_1;
    p.k_degrBV_2=k_degrBV_2;

    p.k_lys=k_lys;
    
    p.randV2_I1=randV2_I1;
    p.randV1_I2=randV1_I2;

    p.k_d_N1=k_d_N1;
    p.k_d_N2=k_d_N2;

    p.K_s=K_s;
    p.Ys_T=Ys_T;
    p.Ys_I1=Ys_I1; 
    p.Ys_I2=Ys_I2; 
    p.Ys_C=Ys_C;
    
    p.threshold=threshold;
    p.rtol=rtol;
    p.Dtau=Dtau;
    p.hmax0=hmax0;

    p.scaling=scaling;
