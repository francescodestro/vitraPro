function [dxdt,dxbind_dt,dxnucl_dt] = inf_model_CS3(t,x,x_bind,x_nucl,...
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
    ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys)

    dxdt=zeros(1,length(x)); % pre-allocate

    b=r_bleed*D;
    
    % species
    T=x(1); 
    V1=x(2);
    V2=x(3); 
    I1=x(startI1:endI1);
    I2=x(startI2:endI2);
    Co=x(startCo:endCo);
    NV=x(endCo+1);
    S=x(endCo+2);
    DIP_all=x(endCo+3);

    VCD=T+sum(I1)+sum(I2)+sum(Co);

    %% pre-allocate all balances 
    dI1dt=zeros(1,length(I1)); 
    dI2dt=dI1dt; 
    dxbind_dt=zeros(1,length(x_bind));
    dxnucl_dt=zeros(1,length(x_nucl));
    dCodt=zeros(1,length(Co)); 
    dNVdt=0;
    dSdt=0;

    %% pre-allocate release scalars
    rel_v1=0;
    rel_v2=0;

    %% uninfected cells: infection and growth
    bind1=k_bind_I1(1)*T*V1; % V1 infection rate counter
    bind2=k_bind_I2(1)*T*V2; % V2 infection rate counter

    growthT=mu_T*T;%*growth_lim;
    dxdt(1)=growthT-k_deathT*T-(bind1+bind2)+Cin*D/scaling-b*x(1); % target cells balance
    dNVdt=dNVdt+k_deathT*T;

    dI1dt(1)=dI1dt(1)+bind1; % infected cells balance
    dI2dt(1)=dI2dt(1)+bind2; % infected cells balance

    dxbind_dt(startB1)=dxbind_dt(startB1)+bind1; % bound virus added to bound virus balance of infected cells of first age bin
    dxbind_dt(startB2)=dxbind_dt(startB2)+bind2; % bound virus added to bound virus balance of infected cells of first age bin

    %% I1, I2: infection, death, and intracellular species balances
    for j = 1:n_bins_inf
        age=age_bins(j);

%         N2_spec=max(0,x_nucl(startB2+j-1)/x(startI2+j-1)); 
%         dI2dt(j)=dI2dt(j)+mu_I2*I2(j)*min(1,60/(60+N2_spec^2)).*(0.5+...
%             0.5*tanh((-VCD+3.7e6/scaling)/0.02));%*min(1,45/(45+age));

%         dI2dt(j)=dI2dt(j)+mu_I2*I2(j);

        dI2dt(j)=dI2dt(j)+mu_I2*I2(j);%*min(1,15/(15+age));
        
        % death dependency on number of DNA copies
        k_deathI1=k_deathT*(0.5+0.5*tanh((tau_death_I1-age)/0.3))+...
            k_death_I1*max(1,log((x_nucl(j))/(I1(j)+1e-50))).*(0.5+...
            0.5*tanh((age-tau_death_I1)/0.3));
        k_deathI2=k_deathT*(0.5+0.5*tanh((tau_death_I2-age)/0.3))+...
            k_death_I2*max(1,log((x_nucl(startB2+j-1))/(I2(j)+1e-50))).*(0.5+...
            0.5*tanh((age-tau_death_I2)/0.3));

        % I1 and I2 death
        deathI1=k_deathI1*I1(j);
        deathI2=k_deathI2*I2(j);
        dI1dt(j)=dI1dt(j)-deathI1;
        dI2dt(j)=dI2dt(j)-deathI2;

        % V1, V2 bound to I1, I2: consumption for cell death and endocytosis
        Ein1=k_int_I1*x_bind(j);
        Ein2=k_int_I2*x_bind(startB2+j-1);
        dxbind_dt(j)=dxbind_dt(j)-Ein1-k_deathI1*x_bind(j);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-Ein2-...
            k_deathI2*x_bind(startB2+j-1); 
        dxnucl_dt(j)=dxnucl_dt(j)+Ein1*eta_I1-k_deathI1*x_nucl(j)-k_d_N1*x_nucl(j);
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)+Ein2*eta_I2-...
            k_deathI2*x_nucl(startB2+j-1)-k_d_N2*x_nucl(startB2+j-1);

        % co-infection: virus of a different type enters infected cell
        kv1=k_bind_I2(j)*V1;
        kv2=k_bind_I1(j)*V2;
        B1_new=I2(j)*kv1; % new virus 1 binding I2
        B2_new=I1(j)*kv2; % new virus 2 binding I1
        bind1=bind1+B1_new; % update counters of bound virus
        bind2=bind2+B2_new; % update counters of bound virus

        dI1dt(j)=dI1dt(j)-B2_new;  % remove new co-infected cells from I1 bal
        dI2dt(j)=dI2dt(j)-B1_new;  % remove new co-infected cells from I2 bal

        dCodt(ind(1+n_bins*(j-1)))=dCodt(ind(1+n_bins*(j-1)))+B2_new; % add new co-infected cells to Co bal
        dCodt(j)=dCodt(j)+B1_new; % add new co-infected cells to Co bal

        v1_toCoB=kv2*x_bind(j); % V1 binded to I1 moving to co-infected cells
        v2_toCoB=kv1*x_bind(startB2+j-1); % V2 binded to I2 moving to co-infected cells
        v1_toCoN=kv2*x_nucl(j); % E1 in I1 moving to co-infected cells
        v2_toCoN=kv1*x_nucl(startB2+j-1); % E2 in I2 moving to co-infected cells

        dxbind_dt(j)=dxbind_dt(j)-v1_toCoB; % move V1 binded to I1 that just got co-infected from B/I1 to B/Co balance
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-v2_toCoB;  % move V2 binded to I2 that just got co-infected from B/I2 to B/Co balance  
        dxnucl_dt(j)=dxnucl_dt(j)-v1_toCoN; % move V1 binded to I1 that just got co-infected from B/I1 to B/Co balance
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)-v2_toCoN;  % move V2 binded to I2 that just got co-infected from B/I2 to B/Co balance  

        dxbind_dt(startB1Co-1+j)=dxbind_dt(startB1Co-1+j)+B1_new; % update V1 binded to new co-infected cells
        dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))+... % update V1 binded to new co-infected cells
            v1_toCoB;
        dxnucl_dt(ind_scB1Co(1+n_bins*(j-1)))=dxnucl_dt(ind_scB1Co(1+n_bins*(j-1)))+... % update E1 to new co-infected cells
            v1_toCoN;

        dxbind_dt(startB2Co-1+j)=dxbind_dt(startB2Co-1+j)+v2_toCoB; % update V2 binded to new co-infected cells
        dxnucl_dt(startB2Co-1+j)=dxnucl_dt(startB2Co-1+j)+v2_toCoN; % update E2 binded to new co-infected cells
        dxbind_dt(ind_scB2Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB2Co(1+n_bins*(j-1)))+B2_new; % update V2 binded to new co-infected cells
            
        % virus of the same type enters I1 and I2
        I1_uptake=k_bind_I1(j)*I1(j)*V1; % calculate virus uptake
        I2_uptake=k_bind_I2(j)*I2(j)*V2;
        dxbind_dt(j)=dxbind_dt(j)+I1_uptake; % add uptake to binded virus balance
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)+I2_uptake;
        bind1=bind1+I1_uptake; % update counters
        bind2=bind2+I2_uptake;

        % nuclear reactions in I1 and I2: BV DNA replication
        dxnucl_dt(j)=dxnucl_dt(j)+k_repl_I1*x_nucl(j)*f_repl_I1(j); % 
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)+k_repl_I2*x_nucl(startB2+j-1)*f_repl_I2(j); % 

        % release from I1
        rel_v1=rel_v1+k_rel_I1*I1(j)*f_rel_I1(j);
        rel_v2=rel_v2+k_rel_I1*I1(j)*f_rel_I1(j)*randV2_I1;

        % release from I2
        rel_v1=rel_v1+k_rel_I2*I2(j)*f_rel_I2(j)*randV1_I2;
        rel_v2=rel_v2+k_rel_I2*I2(j)*f_rel_I2(j);

        % bleeding
        dI1dt(j)=dI1dt(j)-b*I1(j);
        dI2dt(j)=dI2dt(j)-b*I2(j);
        dxbind_dt(startB1+j-1)=dxbind_dt(startB1+j-1)-b*x_bind(startB1+j-1);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-b*x_bind(startB2+j-1);
        dxnucl_dt(startB1+j-1)=dxnucl_dt(startB1+j-1)-b*x_nucl(startB1+j-1);
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)-b*x_nucl(startB2+j-1);

        % viral genome to nonviable cells
        dxnucl_dt(end-1)=dxnucl_dt(end-1)+k_deathI1*x_nucl(j);
        dxnucl_dt(end)=dxnucl_dt(end)+k_deathI2*x_nucl(startB2+j-1);
        dNVdt=dNVdt+k_deathI1*x(startI1+j)+k_deathI2*x(startI2+j-1);
    end
    for j = n_bins_inf+1:n_bins
        age=age_bins(j);

%         N2_spec=max(0,x_nucl(startB2+j-1)/x(startI2+j-1));
%         dI2dt(j)=dI2dt(j)+mu_I2*I2(j)*min(1,15/(15+N2_spec^2));

       dI2dt(j)=dI2dt(j)+mu_I2*I2(j);%*min(1,15/(15+age));%*min(1,60/(60+N2_spec^2)).*(0.5+...
%             0.5*tanh((-VCD+3.7e6)/0.02));%; 
        % death dependency on number of DNA copies
        k_deathI1=k_deathT*(0.5+0.5*tanh((tau_death_I1-age)/0.3))+...
            k_death_I1*max(1,log((x_nucl(j))/(I1(j)+1e-30))).*(0.5+...
            0.5*tanh((age-tau_death_I1)/0.3));
        k_deathI2=k_deathT*(0.5+0.5*tanh((tau_death_I2-age)/0.3))+...
            k_death_I2*max(1,log((x_nucl(startB2+j-1))/(I2(j)+1e-30))).*(0.5+...
            0.5*tanh((age-tau_death_I2)/0.3));


        % I1 and I2 death
        deathI1=k_deathI1*I1(j);
        deathI2=k_deathI2*I2(j);
        dI1dt(j)=dI1dt(j)-deathI1;
        dI2dt(j)=dI2dt(j)-deathI2;


        % V1, V2 bound to I1, I2: consumption for cell death and endocytosis
        Ein1=k_int_I1*x_bind(j);
        Ein2=k_int_I2*x_bind(startB2+j-1);
        dxbind_dt(j)=dxbind_dt(j)-Ein1-k_deathI1*x_bind(j);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-Ein2-...
            k_deathI2*x_bind(startB2+j-1);
        dxnucl_dt(j)=dxnucl_dt(j)+Ein1*eta_I1-k_deathI1*x_nucl(j)-...
            k_d_N1*x_nucl(j);
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)+Ein2*eta_I2-...
            k_deathI2*x_nucl(startB2+j-1)-k_d_N2*x_nucl(startB2+j-1);

        % nuclear reactions in I1 and I2: BV DNA replication
        dxnucl_dt(j)=dxnucl_dt(j)+k_repl_I1*x_nucl(j)*f_repl_I1(j); % 
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)+k_repl_I2*...
            x_nucl(startB2+j-1)*f_repl_I2(j); % 

        % release from I1
        rel_v1=rel_v1+k_rel_I1*I1(j)*f_rel_I1(j);
        rel_v2=rel_v2+k_rel_I1*I1(j)*f_rel_I1(j)*randV2_I1;

        % release from I2
        rel_v1=rel_v1+k_rel_I2*I2(j)*f_rel_I2(j)*randV1_I2;
        rel_v2=rel_v2+k_rel_I2*I2(j)*f_rel_I2(j);

        % bleeding
        dI1dt(j)=dI1dt(j)-b*I1(j);
        dI2dt(j)=dI2dt(j)-b*I2(j);
        dxbind_dt(startB1+j-1)=dxbind_dt(startB1+j-1)-b*x_bind(startB1+j-1);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-b*x_bind(startB2+j-1);
        dxnucl_dt(startB1+j-1)=dxnucl_dt(startB1+j-1)-b*x_nucl(startB1+j-1);
        dxnucl_dt(startB2+j-1)=dxnucl_dt(startB2+j-1)-b*x_nucl(startB2+j-1);

        % viral genome to nonviable cells
        dxnucl_dt(end-1)=dxnucl_dt(end-1)+k_deathI1*x_nucl(startB1+j-1);
        dxnucl_dt(end)=dxnucl_dt(end)+k_deathI2*x_nucl(startB2+j-1);
        dNVdt=dNVdt+k_deathI1*x(startI1+j-1)+k_deathI2*x(startI2+j-1);
    end
    
    % co-infected cells: death, virions binding, endocitosis
    for k=1:length(Co)


        rel_v1=rel_v1+f_rel_Co(k)*k_rel_Co*Co(k)*(x_nucl(startB1Co-1+k))/...
            (x_nucl(startB1Co-1+k)+x_nucl(startB2Co-1+k)+1e-30);
        rel_v2=rel_v2+f_rel_Co(k)*k_rel_Co*Co(k)*(x_nucl(startB2Co-1+k))/...
            (x_nucl(startB1Co-1+k)+x_nucl(startB2Co-1+k)+1e-30);

        k_deathCo=k_deathT*(0.5+0.5*tanh((tau_death_Co-age_co(k))/0.3))+...
            (k_death_Co*max(1,log((f_death_Co1*x_nucl(startB1Co-1+k)+...
            f_death_Co2*x_nucl(startB2Co-1+k))/(Co(k)+1e-30)))).*...
            (0.5+0.5*tanh((age_co(k)-tau_death_Co)/0.3));
        
        deathCo=k_deathCo*Co(k);
        dCodt(k)=dCodt(k)-deathCo;
        
        % V1 binding to co-infected cells
        Co_uptakeV1=k_bind_Co1(k)*Co(k)*V1;
        % V2 binding to co-infected cells
        Co_uptakeV2=k_bind_Co2(k)*Co(k)*V2;

        Ein1=k_int_Co1*x_bind(startB1Co-1+k);
        Ein2=k_int_Co2*x_bind(startB2Co-1+k);

        dxbind_dt(startB1Co-1+k)=dxbind_dt(startB1Co-1+k)+Co_uptakeV1-...
                    k_deathCo*x_bind(startB1Co-1+k)-Ein1;
        dxbind_dt(startB2Co-1+k)=dxbind_dt(startB2Co-1+k)+Co_uptakeV2-...
                    k_deathCo*x_bind(startB2Co-1+k)-Ein2;

        dxnucl_dt(startB1Co-1+k)=dxnucl_dt(startB1Co-1+k)+Ein1*eta_Co1-...
            k_deathCo*x_nucl(startB1Co-1+k)-k_d_N1*x_nucl(startB1Co-1+k);
        dxnucl_dt(startB2Co-1+k)=dxnucl_dt(startB2Co-1+k)+Ein2*eta_Co2-...
            k_deathCo*x_nucl(startB2Co-1+k)-k_d_N2*x_nucl(startB2Co-1+k);
        
        bind1=bind1+Co_uptakeV1;
        bind2=bind2+Co_uptakeV2;
        
        % nuclear reactions in co-infected cells: BV DNA replication
        dxnucl_dt(startB1Co-1+k)=dxnucl_dt(startB1Co-1+k)+k_repl_Co1*...
            x_nucl(startB1Co-1+k)*f_repl_Co1(k); % 
        dxnucl_dt(startB2Co-1+k)=dxnucl_dt(startB2Co-1+k)+k_repl_Co2*...
            x_nucl(startB2Co-1+k)*f_repl_Co2(k); 

        % bleeding
        dCodt(k)=dCodt(k)-b*Co(k);
        dxbind_dt(startB1Co+k-1)=dxbind_dt(startB1Co+k-1)-b*x_bind(startB1Co+k-1);
        dxnucl_dt(startB1Co+k-1)=dxnucl_dt(startB1Co+k-1)-b*x_nucl(startB1Co+k-1);
        dxbind_dt(startB2Co+k-1)=dxbind_dt(startB2Co+k-1)-b*x_bind(startB2Co+k-1);
        dxnucl_dt(startB2Co+k-1)=dxnucl_dt(startB2Co+k-1)-b*x_nucl(startB2Co+k-1);

        % viral genome to nonviable cells
        dxnucl_dt(end-1)=dxnucl_dt(end-1)+k_deathCo*x_nucl(startB1Co-1+k);
        dxnucl_dt(end)=dxnucl_dt(end)+k_deathCo*x_nucl(startB2Co-1+k);
        dNVdt=dNVdt+k_deathCo*x(startCo-1+k);
    end

    %% virions balance
    dxdt(2)=-bind1-k_degrBV_1*V1+rel_v1-D*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-bind2-k_degrBV_2*V2+rel_v2-D*x(3); % BV rep/cap (#/mL) 

    % nonviable cells: outlet, lysis, viral genome degradation
    dNVdt=dNVdt-(k_lys+b)*NV;
    dxnucl_dt(end-1)=dxnucl_dt(end-1)-(k_lys+b+k_d_N1)*x_nucl(end-1);
    dxnucl_dt(end)=dxnucl_dt(end)-(k_lys+b+k_d_N2)*x_nucl(end);

    % compute total DIPs (infective + noninfective)
    dDIPdt=-bind2+rel_v2*300;

    %% put together derivatives vector
    dxdt(4:end)=[dI1dt dI2dt dCodt dNVdt dSdt dDIPdt];
   
end
