% Case study 2
clear

col='k';
lw=4.2;
line=':';

scaling=1e5;

% Simulation parameters
t_final=24*10;

C0=1.5e6;
C0_nv=0;
Din=0.01;
Sin=15*1e3/scaling; % nmol/mL
MOI_1=1;
MOI_2=0;
age_end_inf=20;
age_viab=140;    
S0=15*1e3/scaling; % nmol/mL
Cin=3e6;  
Dt = 1; % h: control interval

%% Control law
Dtau=.02; % mesh grid
tt=0:Dt:t_final;

DD=ones(1,length(tt))*Din;
CCin=ones(1,length(tt))*Cin;
r_bleed_vect=ones(1,length(tt));
    
%% Initialization
p=initialization_parameters(age_viab,age_end_inf,Dtau,Dt,C0,MOI_1,MOI_2,scaling);

n_bins=p.n_bins;
n_bins_inf=p.n_bins_inf;
t=0;
V0_1=C0*MOI_1; % #/mL
sum_h=0;

p.n_bins=1;
x0=[C0 V0_1 0]/scaling;
x0_bind=0;
x0_nucl=0;
age0_bins=0;
    
% Preallocate solution vectors
xx=zeros(length(tt),1e4+2);
xx_bind=zeros(length(tt),1e4);
xx_nucl=xx_bind;
xx(1,1:length(x0))=x0*scaling;
ages_vect=xx_bind;
S_vect=zeros(1,length(tt));
S_vect(1)=S0*scaling;

tic
%% Simulation
for i = 2:length(tt)

    D=DD(i);
    Cin=CCin(i);
    r_bleed=r_bleed_vect(i);

    [t_out,x,x_bind,x_nucl,S,age_bins] = main_MoCh([t t+Dt], x0, x0_bind, ...
        x0_nucl, S0, D, r_bleed, Cin, Sin, p, age0_bins,scaling);

    % initialize following step
    x0=x;
    x0_bind=x_bind;
    x0_nucl=x_nucl;
    age0_bins=age_bins;
    S0=S;
    
    t=t_out(end);

    % save new results
    index=i;
    xx(index,1:length(x))=x*scaling;
    xx_bind(index,1:length(x_bind))=x_bind*scaling;
    xx_nucl(index,1:length(x_nucl))=x_nucl*scaling;
    ages_vect(index,1:length(age_bins))=age_bins;
    S_vect(index)=S*scaling;

end
toc

%% New plots
I1=xx(:,p.startI1:end);
N1=xx_nucl(:,p.startB1:end);

pos=[2079   519     473     245];

figure(1)
hold on
plot(tt/24,sum(I1,2),'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
set(gcf,'Position',pos)
xlabel('Time [d]')
ylabel('Infected [cell/mL]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(2)
hold on
plot(tt/24,xx(:,2),'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22,'ylim',[0 8e8]) %,'xticklabel',[])
set(gcf,'Position',pos)
ylabel('Virion [PFU/mL]')
xlabel('Time [d]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(3)
hold on
plot(tt/24,S_vect/1e3,'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) %,'xticklabel',[])
set(gcf,'Position',pos)
ylabel('Glucose [mM]')
xlabel('Time [d]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(4)
hold on 
plot(tt/24,xx(:,1),'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) %,'ylim',[0 4e6]) %,'xticklabel',[])
set(gcf,'Position',pos)
ylabel('Uninfected [cell/mL]')
xlabel('Time [d]')
box on


%% infection age distribution
t_plot=[48 96 144 240]+1;
pos2=[ 323   454   560   150];
p.age_bins=0:.1:age_viab-Dtau;
age_edges=[p.age_bins age_viab];

for i=1 %:length(t_plot)

    ages_plot=ages_vect(t_plot(i),:);
    I1_moch=I1(t_plot(i),:);

    I1_plot=zeros(1,length(p.age_bins));

    for kk=2:length(I1_plot)
        age_min=p.age_bins(kk-1);
        age_max=p.age_bins(kk);
        age_indices=and(ages_plot>age_min,ages_plot<=age_max);
        I1_plot(kk)=sum(I1_moch(age_indices));
    end
   
    figure(5+i)
    hold on
    plot(p.age_bins,I1_plot,'linewidth',lw,'Color',col,'LineStyle',line)
    
    set(gcf,'Position',pos2);
    box on
end