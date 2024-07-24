% Script for generating Figure S4
close all
clear
clc

addpath('./sub');

Run_CS_S1

par=exp([-1.4438   -2.9241    1.6241]);

Run_CS_S1_lumped(par)

figure(1)
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S3_I.png')

figure(2)
set(gca,'ylim',[0 3e8])
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S3_V.png')

figure(3)
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S3_S.png')

figure(4)
set(gca,'fontsize',30)
set(gcf,'Position',[  4490         172         574         358])
saveas(gcf,'S3_T.png')

legend('This work','Lumped parameter model')

close(figure(6))
