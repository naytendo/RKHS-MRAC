close all; clear; clc;

astar=-10; bstar=+10;
ar=-20; br=20;
ga=10; gb=10; gf=10;
gx=10; gr=10;

x0=5; xr0=0; epsilonc=1.5; epsilone=1;epsilon3m=0.5; epsilon5m=0.5;
alpha0=0; beta0=0;sigma=0.1; N=12;
kx0=0; kr0=0;

%% RKHS Preparation
basiscenter=linspace(-10,10,N)';
basiscenter12=basiscenter;
hype=1/(2*2^2); l=2; sig=2;

% Exponential Kernel
knlfune=@(x,bsc) exp(-hype*(x-bsc).^2);
% Matern 3/2 Kernel
knlfun3m=@(x,bsc) sig^2.*(1+sqrt(3)./l.*abs(x-bsc)).*exp(-sqrt(3).*abs(x-bsc)./l);
% Matern 5/2 Kernel
knlfun5m=@(x,bsc) sig^2.*(1+sqrt(5)/l.*abs(x-bsc)+5/(3*l^2).*abs(x-bsc).^2).*exp(-sqrt(5).*abs(x-bsc)./l);
k_xi_xie=zeros(N);
k_xi_xi3m=zeros(N);
k_xi_xi5m=zeros(N);
for i=1:N
    for j=1:N
        k_xi_xie(i,j)=knlfune(basiscenter(i),basiscenter(j));
    end
end
condnume=cond(k_xi_xie);

inv_k_xi_xie=inv(k_xi_xie);

for i=1:N
    for j=1:N
        k_xi_xi3m(i,j)=knlfun3m(basiscenter(i),basiscenter(j));
    end
end
condnum3m=cond(k_xi_xi3m);

inv_k_xi_xi3m=inv(k_xi_xi3m);

for i=1:N
    for j=1:N
        k_xi_xi5m(i,j)=knlfun5m(basiscenter(i),basiscenter(j));
    end
end
condnum5m=cond(k_xi_xi5m);

inv_k_xi_xi5m=inv(k_xi_xi5m);


%% Run Simulink
open ScalarModel_Deadzone.slx
simresult_Deadzone12=sim('ScalarModel_Deadzone.slx');

open ScalarModel.slx
simresult12=sim('ScalarModel.slx');

%% Post Processing
Theta_hat_E_Dead12=simresult_Deadzone12.Theta_E.Data;
Theta_hat_3m_Dead12=simresult_Deadzone12.Theta_3m.Data;
Theta_hat_5m_Dead12=simresult_Deadzone12.Theta_5m.Data;
AlphaE_Dead12=simresult_Deadzone12.AlphaE.Data;
Alpha3m_Dead12=simresult_Deadzone12.Alpha3m.Data;
Alpha5m_Dead12=simresult_Deadzone12.Alpha5m.Data;
Err_CDead12=simresult_Deadzone12.Err_C.Data;
Err_3mDead12=simresult_Deadzone12.Err_3m.Data;
tsim_Dead12=simresult_Deadzone12.tout;

Theta_hat_E12=simresult12.Theta_E.Data;
Theta_hat_3m12=simresult12.Theta_3m.Data;
Theta_hat_5m12=simresult12.Theta_5m.Data;
AlphaE12=simresult12.AlphaE.Data;
Alpha3m12=simresult12.Alpha3m.Data;
Alpha5m12=simresult12.Alpha5m.Data;
Err_C12=simresult12.Err_C.Data;
Err_3m12=simresult12.Err_3m.Data;
tsim12=simresult12.tout;

%save('Err_Normal','Err_CDead12','Err_3mDead12','Err_C12','Err_3m12','tsim12','tsim_Dead12');
save('Err_Alpha_off','Err_CDead12','Err_3mDead12','Err_C12','Err_3m12','tsim12','tsim_Dead12');

%% Run simulink again with N=24
N=24;
basiscenter=linspace(-10,10,N)';
basiscenter24=basiscenter;
hype=1/(2*2^2); l=2; sig=2;

k_xi_xie=zeros(N);
k_xi_xi3m=zeros(N);
k_xi_xi5m=zeros(N);
for i=1:N
    for j=1:N
        k_xi_xie(i,j)=knlfune(basiscenter(i),basiscenter(j));
    end
end
condnume=cond(k_xi_xie);

inv_k_xi_xie=inv(k_xi_xie);

for i=1:N
    for j=1:N
        k_xi_xi3m(i,j)=knlfun3m(basiscenter(i),basiscenter(j));
    end
end
condnum3m=cond(k_xi_xi3m);

inv_k_xi_xi3m=inv(k_xi_xi3m);

for i=1:N
    for j=1:N
        k_xi_xi5m(i,j)=knlfun5m(basiscenter(i),basiscenter(j));
    end
end
condnum5m=cond(k_xi_xi5m);

inv_k_xi_xi5m=inv(k_xi_xi5m);

simresult_Deadzone24=sim('ScalarModel_Deadzone.slx');
simresult24=sim('ScalarModel.slx');

%% Final Post processing
Theta_hat_E_Dead24=simresult_Deadzone24.Theta_E.Data;
Theta_hat_3m_Dead24=simresult_Deadzone24.Theta_3m.Data;
Theta_hat_5m_Dead24=simresult_Deadzone24.Theta_5m.Data;
AlphaE_Dead24=simresult_Deadzone24.AlphaE.Data;
Alpha3m_Dead24=simresult_Deadzone24.Alpha3m.Data;
Alpha5m_Dead24=simresult_Deadzone24.Alpha5m.Data;
tsim_Dead24=simresult_Deadzone24.tout;

Theta_hat_E24=simresult24.Theta_E.Data;
Theta_hat_3m24=simresult24.Theta_3m.Data;
Theta_hat_5m24=simresult24.Theta_5m.Data;
AlphaE24=simresult24.AlphaE.Data;
Alpha3m24=simresult24.Alpha3m.Data;
Alpha5m24=simresult24.Alpha5m.Data;
tsim24=simresult24.tout;

xplot=[-10:0.1:10];

comp_este_dead12=0; comp_est3m_dead12=0; comp_est5m_dead12=0;
comp_este12=0; comp_est3m12=0; comp_est5m12=0;
for i=1:12
    comp_este_dead12=comp_este_dead12+Theta_hat_E_Dead12(end,i)*knlfune(xplot,basiscenter12(i));
end

for i=1:12
    comp_est3m_dead12=comp_est3m_dead12+Theta_hat_3m_Dead12(end,i)*knlfun3m(xplot,basiscenter12(i));
end

for i=1:12
    comp_est5m_dead12=comp_est5m_dead12+Theta_hat_5m_Dead12(end,i)*knlfun5m(xplot,basiscenter12(i));
end

for i=1:12
    comp_este12=comp_este12+Theta_hat_E12(end,i)*knlfune(xplot,basiscenter12(i));
end

for i=1:12
    comp_est3m12=comp_est3m12+Theta_hat_3m12(end,i)*knlfun3m(xplot,basiscenter12(i));
end

for i=1:12
    comp_est5m12=comp_est5m12+Theta_hat_5m12(end,i)*knlfun5m(xplot,basiscenter12(i));
end

comp_este_dead24=0; comp_est3m_dead24=0; comp_est5m_dead24=0;
comp_este24=0; comp_est3m24=0; comp_est5m24=0;
for i=1:24
    comp_este_dead24=comp_este_dead24+Theta_hat_E_Dead24(end,i)*knlfune(xplot,basiscenter24(i));
end

for i=1:24
    comp_est3m_dead24=comp_est3m_dead24+Theta_hat_3m_Dead24(end,i)*knlfun3m(xplot,basiscenter24(i));
end

for i=1:24
    comp_est5m_dead24=comp_est5m_dead24+Theta_hat_5m_Dead24(end,i)*knlfun5m(xplot,basiscenter24(i));
end

for i=1:24
    comp_este24=comp_este24+Theta_hat_E24(end,i)*knlfune(xplot,basiscenter24(i));
end

for i=1:24
    comp_est3m24=comp_est3m24+Theta_hat_3m24(end,i)*knlfun3m(xplot,basiscenter24(i));
end

for i=1:24
    comp_est5m24=comp_est5m24+Theta_hat_5m24(end,i)*knlfun5m(xplot,basiscenter24(i));
end

%save('Est_Normal','comp_este12','comp_este24','comp_est3m12','comp_est3m_dead12');
save('Est_Alpha_off','comp_este12','comp_este24','comp_est3m12','comp_est3m_dead12');

f=figure(1);
f.Position=[570 50 850 750];
hold on;
plot(xplot,10*tanh(xplot),'k','LineWidth',4);
plot(xplot,comp_este_dead12,'r','LineWidth',2,'Marker','p','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est3m_dead12,'Color','#EDB120','LineWidth',2,'Marker','p','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est5m_dead12,'b','LineWidth',2,'Marker','p','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_este12,'r','LineWidth',2,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est3m12,'Color','#EDB120','LineWidth',2,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_est5m12,'b','LineWidth',2,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',10);
plot(xplot,comp_este_dead24,'r','LineWidth',2,'Marker','s','MarkerIndices',1:10:length(xplot),'Markersize',20);
% plot(xplot,comp_est3m_dead24,'g','LineWidth',1,'Marker','s','MarkerIndices',1:10:length(xplot),'Markersize',15);
% plot(xplot,comp_est5m_dead24,'b','LineWidth',1,'Marker','s','MarkerIndices',1:10:length(xplot),'Markersize',15);
% plot(xplot,comp_este24,'r','LineWidth',1,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',15);
% plot(xplot,comp_est3m24,'g','LineWidth',1,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',15);
% plot(xplot,comp_est5m24,'b','LineWidth',1,'Marker','o','MarkerIndices',1:10:length(xplot),'Markersize',15);

legend('Target Function','Case 1','Case 2','Case 3','Case 4','Case 5',...
    'Case 6','Case 7','Location','southoutside','FontSize',14,'NumColumns',4)

% legend('Target Function','Exponential Kernel Estimation with Deadzone Modification with N=12',...
%     '3/2 Matern Kernel Estimation with Deadzone Modification with N=12',...
%     '5/2 Matern Kernel Estimation with Deadzone Modification with N=12',...
%     'Exponential Kernel Estimation with Sigma Modification with N=12',...
%     '3/2 Matern Kernel Estimation with Sigma Modification with N=12',...
%     '5/2 Matern Kernel Estimation with Sigma Modification with N=12',...
%     'Exponential Kernel Estimation with Deadzone Modification with N=24','Location','southoutside','FontSize',14)%,...
%     '3/2 Matern Kernel Estimation with Deadzone Modification with N=24',...
%     '5/2 Matern Kernel Estimation with Deadzone Modification with N=24',...
%     'Exponential Kernel Estimation with Sigma Modification with N=24',...
%     '3/2 Matern Kernel Estimation with Sigma Modification with N=24',...
%     '5/2 Matern Kernel Estimation with Sigma Modification with N=24','Location','southoutside');
