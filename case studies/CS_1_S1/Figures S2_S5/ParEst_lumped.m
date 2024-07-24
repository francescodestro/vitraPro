% Script for estimating the parameters of the lumped parameter model used
% in Figures S4-S5
clear
clc

addpath('./sub');

load('data_for_lumped.mat')

scaling=1e8;

t_final=24*10;
C0=1.5e6;
C0_nv=0;
S0=15*1e3; % nmol/mL
Cin=3e6;
Din=0.0;
Sin=15*1e3/scaling; % nmol/mL
MOI_1=1;
MOI_2=0;
age_end_inf=20;
age_viab=140;    
Dt = 1; % h: control interval

%% Control law
Dtau=.1; % mesh grid

%% Initialization
p=initialization_parameters(age_viab,age_end_inf,Dtau,Dt,C0,MOI_1,MOI_2,scaling);

%% Parameters lumped model
mu=p.mu_T;
k_bind=p.k_bind0_I1;
k_deathT=p.k_deathT;
k_deathI=p.k_death_I1;
k_rel=p.k_rel_I1;
K_s=p.K_s;
Ys_T=p.Ys_T;
k_degrV=p.k_degrBV_1;

%% Guess parameters
p0=log([k_bind,k_deathI,k_rel]);

%% Estimation
lb=log([.001 0.00001 0.1]); 
ub=log([100 10 50]);

options=optimoptions('fmincon','Algorithm','sqp','UseParallel',true);

p=fmincon(@f_obj,p0,[],[],[],[],lb,ub,[],options,T,V,I,S);

function SSE=f_obj(par,T_exp,V_exp,I_exp,S_exp)

    par=exp(par);

    D=0;
    Cin=0;
    r_bleed=1;
    Sin=0;
    scaling=1e8;
    t_final=24*10;

    t=0;
    C0=1.5e6;
    MOI_1=1;
    V0_1=C0*MOI_1; % #/mL
    S0=15*1e3; % nmol/mL
    Dt = 1; % h: control interval

    x0=[C0 V0_1 0 S0]/scaling;
    
    % Preallocate solution vectors
    tt=0:Dt:t_final;
    xx=zeros(length(tt),4);
    xx(1,:)=x0*scaling;

    mu=0.028;%par(1);
    k_deathT=8e-5;%par(3);
    K_s=1.3e-5;
    Ys_T=1.2e-4;%par(7);
    Ys_I=0;
    k_degrV=0.007;

    k_bind=par(1);
    k_deathI=par(2);
    k_rel=par(3);

    for i = 2:length(tt)
   
        [t_out,x] = ode23s(@inf_model_lumped,[t t+Dt],x0,[],mu,k_bind,k_deathT,...
            k_deathI,k_rel,K_s, Ys_T, Ys_I, k_degrV, Cin, r_bleed, Sin, D, scaling);
    
        % initialize following step
        x0=x(end,:);
        
        t=t_out(end);
    
        % save new results
        index=i;
        xx(index,:)=x(end,:)*scaling;
    
    end

    T=xx(:,1);
    V=xx(:,2);
    I=xx(:,3);
    S=xx(:,4);

    SSE=sum((T-T_exp).^2)/max(T_exp)+sum((V-V_exp).^2)/max(V_exp)+...
        sum((I-I_exp).^2)/max(I_exp)+sum((S-S_exp).^2)/max(S_exp);

end