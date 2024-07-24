% Script for generating Figure S5

close all
clear
clc

addpath('./sub');

Run_CS_S1
Run_CS_S1_moch
Run_CS_S1_FD

figure(1)
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S5_I.png')

figure(2)
set(gca,'ylim',[0 3e8])
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S5_V.png')

figure(3)
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S5_S.png')

figure(4)
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S5_T.png')

figure(6)
set(gca,'ytick',[0 5e4 10e4 15e4])
set(gca,'ylim',[0 15e4],'xlim',[0 30])
set(gca,'fontsize',28)
xlabel('Infection age [hpi]')
ylabel('Infected [cell/mL]')
set(gca,'YAxisLocation','left','linewidth',2,'xtick',0:5:30,...
    'xticklabel',0:5:30,'fontsize',26,'XMinorTick','on')
text(1, 15e4*.9,'1 d post inoculation','FontSize',24)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S5_dist.png')

legend('This work','Method of characteristics','Finite differences')