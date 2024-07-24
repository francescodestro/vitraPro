% Script for generating Figure S3
close all
clear
clc

addpath('./sub');

Run_CS1
Run_CS1_moch
Run_CS1_FD

figure(1)
set(gcf,'Position',[  4490         172         574         358])
set(gca,'fontsize',30)
saveas(gcf,'S4_I.png')

figure(2)
set(gca,'ylim',[0 .8e8])
set(gcf,'Position',[4490    172    574    358])
set(gca,'fontsize',30)
saveas(gcf,'S4_V.png')

figure(3)
set(gcf,'Position',[4490         172         574         358])
set(gca,'fontsize',30)
saveas(gcf,'S4_S.png')

figure(4)
set(gcf,'Position',[4490         172         574         358])
set(gca,'fontsize',30)
saveas(gcf,'S4_T.png')

figure(6)
set(gca,'ytick',[0 5e4 10e4 15e4])
set(gca,'ylim',[0 8e4],'xlim',[20 60])
set(gca,'fontsize',30)
xlabel('Infection age [hpi]')
ylabel('Infected [cell/mL]')
set(gca,'YAxisLocation','left','linewidth',2,'xtick',20:10:60,...
    'xticklabel',20:10:60,'fontsize',26,'XMinorTick','on')
text(21, 8e4*.9,'2 d post inoculation','FontSize',24)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S4_dist.png')

legend('This work','Method of characteristics','Finite differences')