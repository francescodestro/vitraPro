function [dxdt,dxbind_dt,dxnucl_dt] = inf_model(t,x,x_bind,x_nucl,...
    r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, ...
    mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, ...
    f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, ...
    tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, ...
    eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, ...
    k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, ...
    f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, ...
    f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,...
    k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo, ...
    endCo, startB1,startB2, n_bins_inf,age_bins,...
    ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling)

    dxdt=zeros(1,length(x)); % pre-allocate
    b=r_bleed*D;

    % states
    T=x(1);
    V1=x(2); % STV
    I1=x(startI1:endI1);
    NV=x(end-1);
    S=x(end);

    %% pre-allocate all balances 
    dI1dt=zeros(1,length(I1)); % STV-infected
    dxbind_dt=zeros(1,length(x_bind));
    dxnucl_dt=dxbind_dt;
    dNVdt=0;
    dSdt=0;

    %% pre-allocate release scalars
    rel_v1=0;

    %% infection of uninfected cells
    bind1=k_bind_I1(1)*T*V1; % V1 infection rate counter
    growth_lim=max(S/(K_s+S),0);
    growthT=mu_T*T*growth_lim;
    dxdt(1)=growthT-k_deathT*T-(bind1)+Cin*D/scaling-b*x(1); % target cells balance
    dNVdt=dNVdt+k_deathT*T;
    dI1dt(1)=dI1dt(1)+bind1; % infected cells balance
    dxbind_dt(startB1)=dxbind_dt(startB1)+bind1; % bound virus added to bound virus balance of infected cells of first age bin
    
    %% Infected cells growth
    growthI1=mu_I1*I1*growth_lim;

    %% Substrate balance
    dI1dt=dI1dt+growthI1;
%     dSdt=D*(Sin-S)-Ys_T*T-Ys_I1*sum(I1);
    dSdt=D*(Sin-S)-(Ys_T*T+Ys_I1*sum(I1))*S/(1e-2/scaling+S);
    
    %% I1, I2: infection, death, and intracellular species balances
    for j = 1:n_bins_inf
        age=age_bins(j);
        
        % death dependency on number of DNA copies
        k_deathI1=k_deathT*(0.5+0.5*tanh((tau_death_I1-age)/0.3))+...
            k_death_I1*max(1,log((x_nucl(j))/(I1(j)+1e-50))).*(0.5+...
            0.5*tanh((age-tau_death_I1)/0.3));

        % I1 and I2 death
        deathI1=k_deathI1*I1(j);
        dI1dt(j)=dI1dt(j)-deathI1;
        dNVdt=dNVdt+deathI1;

        % V1, V2 binded to I1, I2: consumption for cell death and endocytosis
        Ein1=k_int_I1*x_bind(j);
        dxbind_dt(j)=dxbind_dt(j)-Ein1-k_deathI1*x_bind(j);
        dxnucl_dt(j)=dxnucl_dt(j)+Ein1*eta_I1-(k_deathI1+k_d_N1)*x_nucl(j);

        % virus of the same type enters I1 and I2
        I1_uptake=k_bind_I1(j)*I1(j)*V1; % calculate virus uptake
        dxbind_dt(j)=dxbind_dt(j)+I1_uptake; % add uptake to binded virus balance
        bind1=bind1+I1_uptake; % update counters

        % nuclear reactions in I1 and I2: BV DNA replication
        dxnucl_dt(j)=dxnucl_dt(j)+k_repl_I1*x_nucl(j)*f_repl_I1(j); % 

        % release from I1
        rel_v1=rel_v1+k_rel_I1*I1(j)*f_rel_I1(j);

        % bleeding
        dI1dt(j)=dI1dt(j)-b*I1(j);
        dxbind_dt(startB1+j-1)=dxbind_dt(startB1+j-1)-b*x_bind(startB1+j-1);
        dxnucl_dt(startB1+j-1)=dxnucl_dt(startB1+j-1)-b*x_nucl(startB1+j-1);


    end
    for j = n_bins_inf+1:n_bins
        age=age_bins(j);
        
        % death dependency on number of DNA copies
        k_deathI1=k_deathT*(0.5+0.5*tanh((tau_death_I1-age)/0.3))+...
            k_death_I1*max(1,log((x_nucl(j))/(I1(j)+1e-50))).*(0.5+...
            0.5*tanh((age-tau_death_I1)/0.3));

        % I1 and I2 death
        deathI1=k_deathI1*I1(j);
        dI1dt(j)=dI1dt(j)-deathI1;
        dNVdt=dNVdt+deathI1;

        % V1, V2 binded to I1, I2: consumption for cell death and endocytosis
        Ein1=k_int_I1*x_bind(j);
        dxbind_dt(j)=dxbind_dt(j)-Ein1-k_deathI1*x_bind(j);
        dxnucl_dt(j)=dxnucl_dt(j)+Ein1*eta_I1-(k_deathI1+k_d_N1)*x_nucl(j);

        % nuclear reactions in I1 and I2: BV DNA replication
        dxnucl_dt(j)=dxnucl_dt(j)+k_repl_I1*x_nucl(j)*f_repl_I1(j); % 

        % release from I1
        rel_v1=rel_v1+k_rel_I1*I1(j)*f_rel_I1(j);

        % bleeding
        dI1dt(j)=dI1dt(j)-b*I1(j);
        dxbind_dt(startB1+j-1)=dxbind_dt(startB1+j-1)-b*x_bind(startB1+j-1);
        dxnucl_dt(startB1+j-1)=dxnucl_dt(startB1+j-1)-b*x_nucl(startB1+j-1);
    end
    

    %% virions balance
    dxdt(2)=-bind1-k_degrBV_1*V1+rel_v1-D*x(2); % BV ITR/GOI (#/mL) 

    % outlet nonviable cells
    dNVdt=dNVdt-b*NV;

    %% put together derivatives vector
    dxdt(3:end)=[dI1dt dNVdt dSdt];

end
