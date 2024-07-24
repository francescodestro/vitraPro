% Case study 2
clear

col=[200 20 100]/255;
lw=3;
line='--';

scaling=1e7;

%% Simulation parameters
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
Dt = 1; % h: control interval

%% Control law
Dtau=.1; % mesh grid
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
x0=[C0 V0_1 zeros(1,n_bins) 0 S0]/scaling;
x0_bind=zeros(1,n_bins);
x0_nucl=x0_bind;
sum_h=0;

    
% Preallocate solution vectors
xx=zeros(length(tt),2+n_bins+2);
xx_bind=zeros(length(tt),p.endB1);
xx_nucl=xx_bind;
xx(1,:)=x0*scaling;

tic

%% Simulation
for i = 2:length(tt)

    D=DD(i);
    Cin=CCin(i);
    r_bleed=r_bleed_vect(i);

    [t_out,x,x_bind,x_nucl,sum_h] = main_FD([t t+Dt], x0, x0_bind, ...
        x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling);

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


%% New plots
I1=xx(:,p.startI1:end-2);
N1=xx_nucl(:,p.startB1:end-2);

pos=[2079   519     473     245];

figure(1)
hold on
plot(tt/24,sum(I1,2),'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) 
set(gcf,'Position',pos)
xlabel('Time [d]')
ylabel('Infected [cell/mL]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(2)
hold on
plot(tt/24,xx(:,2),'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22,'ylim',[0 8e8]) 
set(gcf,'Position',pos)
ylabel('Virion [PFU/mL]')
xlabel('Time [d]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(3)
hold on
plot(tt/24,xx(:,end)/1e3,'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) 
set(gcf,'Position',pos)
ylabel('Glucose [mM]')
xlabel('Time [d]')
set(gca,'Units','normalized','OuterPosition',[0 0 1 1])
box on

figure(4)
hold on 
plot(tt/24,xx(:,1),'linewidth',lw,'Color',col,'LineStyle',line)
set(gca,'linewidth',2,'fontsize',22) 
set(gcf,'Position',pos)
ylabel('Uninfected [cell/mL]')
xlabel('Time [d]')
box on

%% infection age distribution
t_plot=[48 96 144 240]+1;
pos2=[323   454   560   150];
for i=1:1
    figure(5+i)
    hold on
    I1_plot=I1(t_plot(i),:);
    plot(p.age_bins,I1_plot,'linewidth',lw,'Color',col,'LineStyle',line)
    set(gca,'YAxisLocation','right','linewidth',2,'xticklabel',[],...
        'fontsize',22,'ylim',[0 1.5e5],'xlim',[0 100],'ytick',[0 5 10]*1e4,...
        'YScale','linear') 
    set(gcf,'Position',pos2);
    box on
end
set(gca,'linewidth',2,'xticklabel',0:20:140,'fontsize',22,'ylim',[0 1.5e5],'xlim',[0 100],'YScale','linear') %,'xticklabel',[])
   