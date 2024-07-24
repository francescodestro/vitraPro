%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%  vitraPro: Simulation of viral transduction and propagation in suspension cultures
%
%       v1.0 April 4, 2024. F. Destro
%
%       Reference:  
%           Destro F. and R.D. Braatz (2024). Efficient simulation of viral
%           transduction and propagation for biomanufacturing. Submitted
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Outputs:
% 
%       Scalar states:
%           --> row i corresponds to values at time tt(i)
%           T:              Uninfected cell concentration time profile - [cell/mL]        
%           V1:             Virion 1 concentration time profile - [PFU/mL]            
%           NV:             Nonviable cell concentration time profile  - [cell/mL]
%           S:              Substrate concentration time profile - [nmol/mL]
%           N1_pc_avg:      Time profile of avg concentration of virus 1 genome in 
%                           infected cells - [vg/cell] (per cell basis)
%
%       States distributed with respect to one infection age: 
%           --> element (i,j) corresponds to cells at infection age equal to
%               p.age_bins(j) at time tt(i)
%           I1:         Time profile of concentration of cells infected by 
%                       virus 1 - [cell/mL] 
%           B1:         Time profile of concentration of virus 1 genome 
%                       bound to infected cells - [vg/mL] (total conc. in system)
%           B1_pc:      Time profile of concentration of virus 1 genome 
%                       bound to infected cells - [vg/mL] (total conc. in system)
%           N1:         Time profile of concentration of virus 1 genome in 
%                       nucleus of infected cells - [vg/mL] (total conc. in system)
%           N1_pc:      Time profile of concentration of virus 1 genome in 
%                       nucleus of infected cells - [vg/mL] (total conc. in system)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Simulation inputs setup
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_final=24*10; % simulation duration [h]
C0_v=1.5e6; % viable cell density at t=0 [cell/mL]
C0_nv=0; % nonviable cell density at t=0 [cell/mL]
S0=15*1e3; % substrate concentration [nmol/mL]

% Inoculation of virions
MOI_1=1; % [PFU/viable cell]

% Inoculation of infected cells
I0_1=0*C0_v/100; % [#/mL] - concentration of inoculated cells infected by virus 1
N_I0_1=5e4*I0_1; % [vg/cell] - viral genome copy number in inoculated cells infected by virus 1
age_I0_1=40; % [hpi] - infection age of inoculated cells infected by virus 

% Discretization parameters
Dtau=.2; % [hpi] mesh spacing - max 2 decimals
age_viab=120;  % [hpi] - maximum infection age in the mesh [hpi] 
               %         select as infection age with low viability and no
               %         recombinant product expression. Must be a multiple
               %         of Dtau and Dt
age_end_inf=10; % [hpi] - infection age beyond which no viral binding occurs.
               %         Must be a multiple of Dtau and Dt. If viral binding occurs 
               %         at all infection ages, set age_end_inf=age_viab     

% Inlet conditions (see Table 1 in Destro and Braatz, 2024)
%   for implementing time-variable inlet conditions, modify the section
%   "Inlet profiles" below
Cin=3e6; % [cell/mL] - viable uninfected cells in feed
Din=0.01; % [1/h] - dilution rate
Sin=15*1e3; % [nmol/mL] - substrate concentration in feed
r_bleed=1; % [–] - bleeding ratio. perfusion: 0 ≤ r < 1; continuous: r = 1

% Inlet profiles 
%   default: constant inlet conditions, defined in previous section
Dt = 1; % [h] control and sampling interval. Must be a multiple of Dtau
tt=0:Dt:t_final;
CCin=ones(1,length(tt))*Cin;
DD=ones(1,length(tt))*Din;
r_bleed_vect=ones(1,length(tt))*r_bleed;
Sin_vect=ones(1,length(tt))*Sin;

% numeric_scheme='RK45';
numeric_scheme='RK23';

generate_figures=1; % 1: generate sample figures; 0: do not generate sample figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation (don't modify)
scaling=5e7;
p=def_parameters(age_viab,age_end_inf,Dtau,scaling);

n_bins=p.n_bins;
n_bins_inf=p.n_bins_inf;
t=0;
V0_1=C0_v*MOI_1; % PFU/mL
x0=[C0_v V0_1 zeros(1,n_bins) C0_nv S0]/scaling;
x0_bind=zeros(1,p.endB1);
x0_nucl=[x0_bind 0];
sum_h=0;
    
% Preallocate solution vectors
xx=zeros(length(tt),length(x0));
xx_bind=zeros(length(tt),p.endB1);
xx_nucl=zeros(length(tt),p.endB1+1);
xx(1,:)=x0*scaling;
    
Sin_vect=Sin_vect/scaling;

% Inoculation of infected cells
bin_inoculation_1=round(age_I0_1/Dtau);
x0(p.startI1+bin_inoculation_1)=I0_1/scaling;
x0_nucl(p.startB1+bin_inoculation_1)=N_I0_1/scaling;

tic
for i = 2:length(tt)

    D=DD(i);
    Cin=CCin(i);
    r_bleed=r_bleed_vect(i);
    Sin=Sin_vect(i);

    [t_out,x,x_bind,x_nucl,sum_h] = main([t t+Dt], x0, x0_bind, ...
        x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling, numeric_scheme);

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
toc

%% Simulation outputs: detailed description in function header
T=xx(:,1); % [cell/mL] - scalar
V1=xx(:,2); % [PFU/mL] - scalar
NV=xx(:,end-1); % [cell/mL] - scalar
S=xx(:,end); % [nmol/mL] - scalar

I1=xx(:,p.startI1:p.endI1); % [cell/mL] - distribution
B1=xx_bind(:,p.startB1:p.endB1); % [vg/mL]  - distribution
B1_pc=B1./I1; % [vg/cell] - distribution
N1=xx_nucl(:,p.startB1:p.endB1); % [vg/mL]  - distribution
N1_pc=N1./I1; % [vg/cell] - distribution

N1_pc_avg=sum(N1,2)./sum(I1,2); % [vg/cell] - scalar

%% Plots
if generate_figures
    sample_figures;
end