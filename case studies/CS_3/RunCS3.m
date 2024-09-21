% Run case study 3 and generate Figure 8
%   run figures_SI for generating Figures S6-S8

clc, clear,% close all

% Load experimental data from:
%   Rüdiger, D., Pelz, L., Hein, M.D., Kupke, S.Y. and Reichl, U., 2021. 
%   Multiscale model of defective interfering particle replication for 
%   influenza A virus infection in animal cell culture. PLoS Computational 
%   Biology, 17(9), p.e1009357.
load influenza_data.mat 

%% Simulation
exps=1:12;

cols=[{[189 148 196]/255} {[163 202 235]/255} {[51 74 159]/255} {[173 76 157]/255} {[0 169 192]/255}];
lines=[{'-.'},{'-'},{':'}];
lws=[2 2 3.5];
xlims=[35 35 64 40 35 35 35 35 35 35 35 35];

numeric_scheme='RK23';

tic
parfor i=exps

    VCC=data{i}.VCC;
    VCC=VCC(~isnan(VCC));
    C0_v=mean(VCC(1:2))*(1-mean(data{i}.T_nonv(1:3))/100);
    C0_nv=mean(VCC(1:2))*mean(data{i}.T_nonv(1:3))/100;

    [t,T,I1,I2,Co,STV_PFU,DNA1_v,DNA1_nv,DNA_dip,DIP_inf,DIP_all,NV]=...
        run_exp_CS3(xlims(i)+.75,C0_v,C0_nv,...
         data{i}.MOI_stv, data{i}.MOI_dip, numeric_scheme);

    t=t-0.75; % inoculation at time -0.75 h -> extra simulation time 

    t_sim{i}=t;
    T_sim{i}=T; 
    I1_sim{i}=I1;
    I2_sim{i}=I2;
    Co_sim{i}=Co;
    STV_PFU_sim{i}=STV_PFU;
    DNA1_v_sim{i}=DNA1_v;
    DNA1_nv_sim{i}=DNA1_nv;
    DIP_sim{i}=DIP_all;
    NV_sim{i}=NV;
    DNA2_sim{i}=DNA_dip;
end
toc

save('simulation_results.mat')

%% Generate Figure 8
for i=1:12
    figure(i)
    hold on
    plot(data{i}.t,data{i}.STV_PFU,'o','MarkerSize',8,'Linewidth',lws(1),'Color',cols{4},...
        'LineStyle','none')
    plot(t_sim{i},STV_PFU_sim{i},'Linewidth',lws(1),'Color',cols{4},...
        'LineStyle',lines{1})
    plot(data{i}.t,data{i}.DIP,'v','MarkerSize',8,'Linewidth',lws(1),'Color',cols{3},...
         'LineStyle','none')
    plot(t_sim{i},DIP_sim{i},'Linewidth',lws(1),'Color',cols{3},...
        'LineStyle',lines{2})
    
    set(gca,'LineWidth',1.5,'fontsize',22,'ylim',[10 10e11],'YScale','log',...
        'ytick',[1e1 1e4 1e8 1e12])
    set(gca,'YMinorTick','on')
    box on
    set(gca,'xlim',[0 xlims(i)],'xtick',0:10:100);
    
    set(gcf,'Position',[ 871   272   357   264])
    if i==1 || i==5 || i==9
         ylabel([{'Viral titer'}])
    end
    if i==9 || i==10 || i==11 || i==12
         xlabel('Time [hpi]')
    end
end

x_text_MOI=[1, -.9, 1.9, 0.375, 3, 1.1, 3.5, 2.6, 2, 0.55, 2.8, 1.9];
for i= [1 2 3 4 6 8 10 12]
    figure(i)
    text(0.7143*xlims(i),50,'Training','FontSize',18)
    text(x_text_MOI(i),5e12,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end
for i= [5 7 9 11]
    figure(i)
    text(0.66*xlims(i),50,'Validation','FontSize',18)
    text(x_text_MOI(i),5e12,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')

end
figure(9)
set(gcf,'Position',[  475   250   386   299])
for j=[1 5]
    figure(j)
    set(gcf,'Position',[  629   237   389   264])
end

for j=10:12
    figure(j)
    set(gcf,'Position',[  454   220   357   297])
end

% for i=1:12
%     figure(i)
%     saveas(gcf,['CS3_titer_' num2str(i) '.png'])
% end