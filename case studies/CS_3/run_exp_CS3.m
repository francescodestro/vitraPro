function [tt,T,I1,I2,Co,STV_PFU,DNA_stv_v,DNA_stv_nv,DNA_dip,DIP,DIP_all,NV,p]=...
    run_exp_CS3(t_final,C0_v,C0_nv,...
        MOI_1,MOI_2, numeric_scheme)

    age_viab=round(t_final);
    age_end_inf=age_viab-1;

    Dtau=.25; 
    Dt = .5; 
    scaling=1e7;
    S0=0; 
    
    Cin=0;
    Din=0;
    Sin=0;
    r_bleed=0;

    %% Initialization
    p=def_parameters(age_viab,age_end_inf,Dtau,scaling);

    n_bins=p.n_bins;
    n_bins_inf=p.n_bins_inf;
    t=0;
    V0_1=2.2e6*MOI_1; 
    V0_2=2.2e6*MOI_2; 
    x0=[C0_v V0_1 V0_2 zeros(1,n_bins*2+n_bins^2-(n_bins-n_bins_inf)*(n_bins-...
        n_bins_inf+1)) C0_nv S0 V0_2]/scaling;
    
    x0_bind=zeros(1,p.endB2Co);
    x0_nucl=[x0_bind 0 0];
    
    sum_h=0;
        
    % Preallocate solution vectors
    tt=0:Dt:t_final+1;
    xx=zeros(length(tt),length(x0));
    xx_bind=zeros(length(tt),p.endB2Co);
    xx_nucl=zeros(length(tt),p.endB2Co+2);
    xx(1,:)=x0*scaling;
        
    CCin=ones(1,length(tt))*Cin;
    DD=ones(1,length(tt))*Din;
    r_bleed_vect=ones(1,length(tt))*r_bleed;
    
    tic
    
    %% Simulation
    for i = 2:length(tt)
    
        D=DD(i);
        Cin=CCin(i);
        r_bleed=r_bleed_vect(i);

    switch numeric_scheme

        case 'RK23'
    
        [t_out,x,x_bind,x_nucl,sum_h] = main_CS3_RK23([t t+Dt], x0, x0_bind, ...
            x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling);

        case 'RK45'

        [t_out,x,x_bind,x_nucl,sum_h] = main_CS3_RK45([t t+Dt], x0, x0_bind, ...
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

    %% Results
    T=xx(:,1);
    I1=sum(xx(:,p.startI1:p.endI1),2);
    I2=sum(xx(:,p.startI2:p.endI2),2);
    Co=sum(xx(:,p.startCo:p.endCo),2);
    
    DIP=xx(:,3);
    STV_PFU=xx(:,2);

    NV=xx(:,end-2);

    DNA1_Co=xx_nucl(:,p.startB1Co:p.endB1Co);
    DNA1_I1=xx_nucl(:,p.startB1:p.endB1);
    DNA1_NV=xx_nucl(:,end-1);
    DNA2_Co=xx_nucl(:,p.startB2Co:p.endB2Co);
    DNA2_I2=xx_nucl(:,p.startB2:p.endB2);
    DNA2_NV=xx_nucl(:,end);

    DNA_stv_v=sum(DNA1_I1,2)+sum(DNA1_Co,2);
    DNA_stv_nv=DNA1_NV;
    DNA_dip=sum(DNA2_I2,2)+sum(DNA2_Co,2)+DNA2_NV;

    DIP_all=xx(:,end);
end