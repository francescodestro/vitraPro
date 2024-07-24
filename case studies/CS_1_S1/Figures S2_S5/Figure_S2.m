% Script for generating Figure S2
close all
clear
clc

addpath('./sub');

Run_CS1

par=exp([-1.4438   -2.9241    1.6241]);

Run_CS1_lumped(par)

figure(1)
set(gcf,'Position',[  4490         172         574         358])
set(gca,'fontsize',30)
saveas(gcf,'S2_I.png')

figure(2)
set(gca,'ylim',[0 3e8])
set(gcf,'Position',[4490    172    574    358])
set(gca,'fontsize',30)
saveas(gcf,'S2_V.png')

figure(3)
set(gcf,'Position',[4490         172         574         358])
set(gca,'fontsize',30)
saveas(gcf,'S2_S.png')

figure(4)
set(gcf,'Position',[4490         172         574         358])
set(gca,'fontsize',30)
saveas(gcf,'S2_T.png')

close(figure(6))