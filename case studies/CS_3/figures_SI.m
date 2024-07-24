% Load simulation results for case study 3
clc, clear, close all

load simulation_results.mat

%% Generate Figure S6 (vRNA)
close all
for i=1:12
    col=cols{2};
    figure(i)

    DNA1_sim=DNA1_v_sim{i}+DNA1_nv_sim{i};
    hold on
    plot(data{i}.t,data{i}.FL_vRNA,'o','MarkerSize',8,'Linewidth',lws(1),'Color',cols{4},...
        'LineStyle','none')

    plot(t_sim{i},DNA1_sim./(I1_sim{i}+I2_sim{i}+Co_sim{i}+T_sim{i}+NV_sim{i}),...
        'Linewidth',lws(1),'Color',cols{4},'LineStyle',lines{1})
    set(gca,'LineWidth',1.5,'fontsize',22)
    box on
    set(gca,'xlim',[0 xlims(i)]);
    set(gca,'YScale','log','ylim',[1e-2 5e6],'ytick',[1e-2 100 1e6])
    set(gcf,'Position',[ 871   272   357   264])
    if i==1 || i==5 || i==9
         ylabel([{'STV copy number'},{'[vRNA/cell]'}])
    end
    if i==9 || i==10 || i==11 || i==12
         xlabel('Time [hpi]')
    end
    set(gca,'YMinorTick','on')
    plot(data{i}.t,data{i}.DI_vRNA,'v','MarkerSize',8,'Linewidth',lws(1),'Color',cols{3},...
         'LineStyle','none')
    hold on
    plot(t_sim{i},DNA2_sim{i}./(I1_sim{i}+I2_sim{i}+Co_sim{i}+T_sim{i}+NV_sim{i}),...
        'Linewidth',lws(1),'Color',cols{3},'LineStyle',lines{2})
end
x_text_MOI=[1, -.9, 1.9, 0.375, 3, 1.1, 3.5, 2.6, 2, 0.55, 2.8, 1.9];
for i= [1 2 3 4 6 8 10 12]
    figure(i)
    text(0.7143*xlims(i),.035,'Training','FontSize',18)
    text(x_text_MOI(i),1.5e7,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end
for i= [5 7 9 11]
    figure(i)
    text(0.66*xlims(i),.03,'Validation','FontSize',18)
    text(x_text_MOI(i),1.5e7,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end

figure(9)
set(gcf,'Position',[ 817   327   420   301])

for j=[1 5]
    figure(j)
    set(gcf,'Position',[   380   545   418   264])
end

for j=10:12
    figure(j)
    set(gcf,'Position',[  454   220   357   297])
end

for i=1:12
    figure(i)
    saveas(gcf,['CS3_vRNA_' num2str(i) '.png'])
end
%% Generate Figure S7 (VCD)
close all

for i=1:12
    figure
        t=t_sim{i};
    plot(data{i}.t,data{i}.VCC,'o','MarkerSize',8,'Linewidth',lws(1),'Color',cols{5},...
        'LineStyle','none')
    hold on
    plot(t_sim{i},T_sim{i}+NV_sim{i}+I1_sim{i}+I2_sim{i}+Co_sim{i},'Linewidth',lws(1),'Color',cols{5})
    set(gca,'LineWidth',1.5,'fontsize',22,'ylim',[1e4 1e7],'yscale','log')
    set(gca,'xlim',[0 xlims(i)]);
    box on
    hold on
    set(gcf,'Position',[ 871   272   357   264])
    set(gca,'YMinorTick','on')
    if i==1 || i==5 || i==9
        ylabel([{'Viable cell density'},{'[cell/mL]'}])
    end
    if i==9 || i==10 || i==11 || i==12
        xlabel('Time [hpi]')
    end
end
x_text_MOI=[1, -.9, 1.9, 0.375, 3, 1.1, 3.5, 2.6, 2, 0.55, 2.8, 1.9];
for i= [1 2 3 4 6 8 10 12]
    figure(i)
    text(0.7143*xlims(i),1.5e4,'Training','FontSize',18)
    text(x_text_MOI(i),1e7*1.5,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end
for i= [5 7 9 11]
    figure(i)
    text(0.66*xlims(i),1.5e4,'Validation','FontSize',18)
    text(x_text_MOI(i),1e7*1.5,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end
figure(9)
set(gcf,'Position',[ 817   327   420   301])
for j=[1 5]
    figure(j)
    set(gcf,'Position',[   380   545   418   264])
end
for j=10:12
    figure(j)
    set(gcf,'Position',[  454   220   357   297])
end
for i=1:12
    figure(i)
    saveas(gcf,['CS3_VCD_' num2str(i) '.png'])
end

%% Generate Figure S8 (infected cells)
close all
for i=1:12
    col=cols{2};
    figure
    plot(data{i}.t,data{i}.VCC.*data{i}.I/100,'o','MarkerSize',8,'Linewidth',lws(1),'Color',cols{5},...
        'LineStyle','none')
    hold on
    plot(t_sim{i},I1_sim{i}+I2_sim{i}+Co_sim{i},'Linewidth',lws(1),'Color',cols{5})
    set(gca,'LineWidth',1.5,'fontsize',22,'ylim',[1e1 1e7],'yscale','log',...
        'ytick',[1e1 1e4 1e7])
    set(gca,'xlim',[0 xlims(i)]);
    box on
    hold on
    set(gcf,'Position',[ 871   272   357   264])
    set(gca,'YMinorTick','on')
    if i==1 || i==5 || i==9
         ylabel([{'Infected cells'},{'[cell/mL]'}])
    end
    if i==9 || i==10 || i==11 || i==12
         xlabel('Time [hpi]')
    end
end
x_text_MOI=[1, -.9, 1.9, 0.375, 3, 1.1, 3.5, 2.6, 2, 0.55, 2.8, 1.9];
for i= [1 2 3 4 6 8 10 12]
    figure(i)
    text(0.66*xlims(i),22,'Validation','FontSize',18)
    text(x_text_MOI(i),1e7*2.2,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end
for i= [5 7 9 11]
    figure(i)
    text(0.66*xlims(i),22,'Validation','FontSize',18)
    text(x_text_MOI(i),1e7*2.2,['MOI STV = ' num2str(data{i}.MOI_stv) ', MOI DIP = ' num2str(data{i}.MOI_dip)],'FontSize',18, 'fontweight', 'bold')
end
figure(9)
set(gcf,'Position',[ 817   327   420   301])
for j=[1 5]
    figure(j)
    set(gcf,'Position',[   380   545   418   264])
end
for j=10:12
    figure(j)
    set(gcf,'Position',[  454   220   357   297])
end
for i=1:12
    figure(i)
    saveas(gcf,['CS3_I_' num2str(i) '.png'])
end