clear;
close all;
clc;

xplot=[-10:0.1:10];
load Err_Normal.mat
Err_C=Err_C12;
Err_CDead=Err_CDead12;
Err_3m_Normal=Err_3m12;
Err_3m_Normal_Dead=Err_3mDead12;
tsim_Normal=tsim12;
tsim_Dead_Normal=tsim_Dead12;
load Err_Alpha_off.mat
Err_3m_Alpha_off=Err_3m12;
Err_3m_Alpha_off_Dead=Err_3mDead12;
tsim_Alpha_off=tsim12;
tsim_Dead_Alpha_off=tsim_Dead12;

f1=figure(1);
f1.Position=[200 200 800 400];
hold on
plot(tsim_Dead_Alpha_off,Err_3m_Alpha_off_Dead,'-.','Color','#0072BD','Linewidth',2);
plot(tsim_Alpha_off,Err_3m_Alpha_off,':','Color','#D95319','Linewidth',2);
plot(tsim_Normal,Err_3m_Normal,'Color','#EDB120','Linewidth',2);
grid on

legend('3/2 Matern Dead-Zone Alpha Off','3/2 Matern Sigma Alpha Off','3/2 Matern Sigma')

load Est_Normal.mat
comp_este12_Normal=comp_este12;
comp_este24_Normal=comp_este24;
comp_est3m12_Normal=comp_est3m12;
comp_est3m_dead12_Normal=comp_est3m_dead12;
load Est_Alpha_off.mat
comp_est3m12_Alpha_off=comp_est3m12;
comp_est3m_dead12_Alpha_off=comp_est3m_dead12;

f2=figure(2);
f2.Position=[1000 200 600 500];
hold on;
plot(xplot,10*tanh(xplot),'k','LineWidth',4);
plot(xplot,comp_este12_Normal,'r','LineWidth',2,'Marker','>','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_este24_Normal,'r','LineWidth',2,'Marker','s','MarkerIndices',1:10:length(xplot),'Markersize',15);
plot(xplot,comp_est3m12_Normal,'b','LineWidth',2,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est3m_dead12_Normal,'b','LineWidth',2,'Marker','d','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est3m12_Alpha_off,'b','LineWidth',2,'Marker','+','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est3m_dead12_Alpha_off,'b','LineWidth',2,'Marker','p','MarkerIndices',1:10:length(xplot),'Markersize',10);
legend('Target Function','Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Location','southoutside','FontSize',14,'NumColumns',4)

