close all; clear; clc;

astar=-10; bstar=+10;
ar=-20; br=20;
ga=10; gb=10; gf=10;
gx=10; gr=10;

x0=0; xr0=0; epsilonc=1;epsilone=1;epsilon3m=0.5; epsilon5m=0.5;
alpha0=0; beta0=0;sigma=0.1; N=12;
kx0=0; kr0=0;

%% RKHS Preparation
basiscenter=linspace(-10,10,N)';
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

Thetastare=inv_k_xi_xie*(10*tanh(basiscenter));

for i=1:N
    for j=1:N
        k_xi_xi3m(i,j)=knlfun3m(basiscenter(i),basiscenter(j));
    end
end
condnum3m=cond(k_xi_xi3m);

inv_k_xi_xi3m=inv(k_xi_xi3m);

Thetastar3m=inv_k_xi_xi3m*(10*tanh(basiscenter));

for i=1:N
    for j=1:N
        k_xi_xi5m(i,j)=knlfun5m(basiscenter(i),basiscenter(j));
    end
end
condnum5m=cond(k_xi_xi5m);

inv_k_xi_xi5m=inv(k_xi_xi5m);

Thetastar5m=inv_k_xi_xi5m*(10*tanh(basiscenter));

%% Preliminary results
xplot=[-10:0.1:10];
complete_estimatione=0;
for i=1:N
    complete_estimatione=complete_estimatione+Thetastare(i)*knlfune(xplot,basiscenter(i));
end

complete_estimation3m=0;
for i=1:N
    complete_estimation3m=complete_estimation3m+Thetastar3m(i)*knlfun3m(xplot,basiscenter(i));
end

complete_estimation5m=0;
for i=1:N
    complete_estimation5m=complete_estimation5m+Thetastar5m(i)*knlfun5m(xplot,basiscenter(i));
end

figure(1)
hold on
plot(xplot,10*tanh(xplot),'k','Linewidth',2);
plot(xplot,complete_estimatione,'Linewidth',2);
plot(xplot,complete_estimation3m,'Linewidth',2);
plot(xplot,complete_estimation5m,'Linewidth',2);
legend('Real Funtion','Estimated Function with Exponentional Kernel',...
    'Estimated Function with 3/2 Matern Kernel','Estimated Function with 5/2 Matern Kernel')

figure (2)
hold on
for i=1:N
    plot(xplot,knlfune(xplot,basiscenter(i)),'Linewidth',2);
end
figure (3)
hold on
for i=1:N
    plot(xplot,knlfun3m(xplot,basiscenter(i)),'Linewidth',2);
end
figure (4)
hold on
for i=1:N
    plot(xplot,knlfun5m(xplot,basiscenter(i)),'Linewidth',2);
end

%% Run Simulink
open ScalarModel_Deadzone.slx
simresult=sim('ScalarModel_Deadzone.slx');

%% Post Processing
Theta_hat_E=simresult.Theta_E.Data;
Theta_hat_3m=simresult.Theta_3m.Data;
Theta_hat_5m=simresult.Theta_5m.Data;
tsim=simresult.tout;

ref_surf=zeros(length(tsim),length(xplot));
for i=1:length(tsim)
    ref_surf(i,:)=10*tanh(xplot);
end

hat_estimatione=zeros(length(tsim),length(xplot));
for i=1:length(tsim)
    for j=1:N
        hat_estimatione(i,:)=hat_estimatione(i,:)+Theta_hat_E(i,j)*knlfune(xplot,basiscenter(j));
    end
end

hat_estimation3m=zeros(length(tsim),length(xplot));
for i=1:length(tsim)
    for j=1:N
        hat_estimation3m(i,:)=hat_estimation3m(i,:)+Theta_hat_3m(i,j)*knlfun3m(xplot,basiscenter(j));
    end
end

hat_estimation5m=zeros(length(tsim),length(xplot));
for i=1:length(tsim)
    for j=1:N
        hat_estimation5m(i,:)=hat_estimation5m(i,:)+Theta_hat_5m(i,j)*knlfun5m(xplot,basiscenter(j));
    end
end

figure(5)
hold on
s2=surface(xplot,tsim,ref_surf,'FaceAlpha',0.2,'FaceColor','r');
s1=surface(xplot,tsim,hat_estimatione);
view(3);
s1.EdgeColor='none';
s2.EdgeColor='none';
ylabel('Time(s)');
xlabel('x Value');
zlabel('Function Value from Exponential Kernel')
legend('Unknown Function','Function Estimation')

figure(6)
hold on
s2=surface(xplot,tsim,ref_surf,'FaceAlpha',0.2,'FaceColor','r');
s1=surface(xplot,tsim,hat_estimation3m);
view(3);
s1.EdgeColor='none';
s2.EdgeColor='none';
ylabel('Time(s)');
xlabel('x Value');
zlabel('Function Value from 3/2 Matern Kernel')
legend('Unknown Function','Function Estimation')

figure(7)
hold on
s2=surface(xplot,tsim,ref_surf,'FaceAlpha',0.2,'FaceColor','r');
s1=surface(xplot,tsim,hat_estimation5m);
view(3);
s1.EdgeColor='none';
s2.EdgeColor='none';
ylabel('Time(s)');
xlabel('x Value');
zlabel('Function Value from 5/2 Matern Kernel')
legend('Unknown Function','Function Estimation')