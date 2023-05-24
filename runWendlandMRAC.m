% clear; close all;clc
delta = 0;

ext = 0.1;
A=-4; B=0.8;
Ar=-5; Br=1;
ga = 10; 
gb = 100; 
gf = 1;
fMag = 1;
SNU = 10; %Signal to Uncertainty
OmegaLim = 5;


x0=OmegaLim+5; xr0=OmegaLim+3;
alpha0=0; beta0=0;sigma=0.1;
kx0=0; kr0=0;


% type = 'exp';
type = 'g';
Dim = 1;
scale=0.5;

s = 1.5;
nu = 1/2*(s+Dim/2);
par = 1;

%% Target Function generation
sTarget = 3;
targetNu = 1/2*(sTarget+Dim/2);
Ip = 8;
targetScale = 0.2;
interps = linspace(-OmegaLim,OmegaLim,Ip)';

targetType = 'ms';


% 
dmat = distsqh(interps,interps);

intmat = kermat(interps,interps,targetType,targetNu-Dim/2,targetScale);
ySamps = fMag*tanh(interps);

targetCoefs = intmat\ySamps;

tol = 0;
maxN = 128;
h = 1e-3;
[gCenters,error,powBound,norm_fs,condNum,condNew] = greedyCenters(tol,h,OmegaLim,ext,maxN,type,par,scale,interps,targetType,targetNu,targetScale,targetCoefs,Dim);


% figure(1)
% hold on
% x = -20:0.1:20;
% ySmooth = 10*tanh(x);
% yAct = zeros(length(x),1);
% 
% for ii = 1:length(x)
%     targetKnlV = zeros(Ip,1);
%     for kk = 1:Ip
%         targetKnlV(kk) = kernel('matern',interps(kk),x(ii),targetHyper,targetNu-Dim/2);
%         yAct(ii) = targetCoefs'*targetKnlV;
%     end
% end

% smoothLine = plot(x,ySmooth,'-.','linewidth',1.5);
% plot(x,yAct,'linewidth',1.5)
% plot(interps,ySamps,'ko','linewidth',1.5)
% legend('10tanh(x)','f_s(x)','interpolating points','location','southeast')
% xlabel('x')
% set(gca,'Fontsize',20)
% Ns = [4+1,8+1,16+1,32+1];
% Ns = [4+1,8+1,16+1,32+1,64+1];

runTimes = 50;


% [offlineError,PowerBound,condNum,norm_fs] = getOfflineError(Ns,type,par,scale,OmegaLim,targetNu,targetType,targetScale,targetCoefs,interps,Dim,ext);

model = 'WendlandMRAC';
open(model)

for tt = 1:length(runTimes)

ssErrorWendland = zeros(length(gCenters),1);
approxError = zeros(length(gCenters),1);
condNum = zeros(length(gCenters),1);
condNumKL = zeros(length(gCenters),1);
cs = getActiveConfigSet(model);
csNew = copy(cs);
set_param(csNew,"StopTime",num2str(runTimes(tt)))
StartMeasureTime = runTimes(tt)-1;





for nn = 1:10:length(gCenters)
% N = Ns(nn);
N = nn;
epsilon = powBound(nn)*norm_fs+delta;

%% RKHS Preparation
% centers=linspace(-OmegaLim-ext,OmegaLim+ext,N)';
centers=gCenters(1:nn);
K = kermat(centers,centers,type,par,scale);
condNum(nn)=cond(K);





%% Run Simulink

simresult=sim(model,csNew);

timeSeriesSSError = simresult.ssErrorWendland;
ssErrorWendland(nn) = max(timeSeriesSSError.Data);

timeSeriesApproxError = simresult.approxError;
approxError(nn) = max(timeSeriesApproxError.Data);

% x = -8.5:0.01:8.5;
% kv = zeros(length(centers),1);
% Pow = zeros(length(x),1);
% 
% for ii = 1:length(x)
%     for kk = 1:length(centers)
%         kv(kk) = kernel(type,centers(kk),x(ii),hyper,nu);
%     end
%     Pow(ii) = powerFunction(type,x,kv,K,hyper);
% 
% end
% PowerBound(nn) = max(Pow);
% figure(5)
% hold on
% plot(x,Pow,'linewidth',1.5)

% figure()
% trackingError = simresult.trackingError;
% plot(trackingError.time,trackingError.Data,'linewidth',1.5)
% hold on
% plot(trackingError.time,epsilon*ones(length(trackingError.time),1),'k--','linewidth',1.5)
% 
% plot(trackingError.time,-epsilon*ones(length(trackingError.time),1),'k--','linewidth',1.5)
% legend('$\tilde{x}_N(t)$','$\epsilon +\delta$ bound','interpreter','latex')
% xlabel('time (s)')
% ylabel('$\tilde{x}_N$','interpreter','latex')
% 
% set(gca,'fontsize',20)
end

%%
% figure(2)
% hold on
% loglog(fillDist,ssErrorWendland,'-o','linewidth',1.5)
% % loglog(fillDist,0.5*fillDist.^0.5,'k-.','linewidth',1.5)
% xlabel('$h_{\Xi_n,M}$','interpreter','latex')
% ylabel('Steady State error','interpreter','latex')
% set(gca,'fontsize',20)
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% 
% figure(3)
% 
% hold on
% loglog(Ns,condNum,'-o','linewidth',1.5)
% xlabel('$N$','interpreter','latex')
% ylabel('Condition Number','interpreter','latex')
% set(gca,'fontsize',20)
% set(gca,'XScale','log')
% set(gca,'YScale','log')

fSS = figure(4);
hold on
lineSS = loglog(1:10:length(gCenters),ssErrorWendland(1:10:length(gCenters)),'-o','linewidth',2.5,'color',[0.9290 0.6940 0.1250]);
colorSS = get(lineSS,'color');

% loglog(Ns,PowerBound,'k-.','linewidth',1.5);


% fApprox = figure(12);
% hold on
% loglog(Ns,approxError,'-o','linewidth',1.5)


% loglog(Ns,PowerBound,'k-.','linewidth',1.5);


end
%%
figure(fSS)
% loglog(Ns,PowerBound+delta,'k--','linewidth',1.5);

% loglog(Ns,delta*ones(length(Ns),1),'k-.','linewidth',1.5)
xlabel('$N$','interpreter','latex')
ylabel('$\frac{\|\tilde{x}_N\|}{\|B^TP\|}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')
% legend('$\frac{\|\tilde{x}_N\|}{\|B^TP\|}$ @ SS','$\epsilon$','$\delta$','interpreter','latex')
% 
% figure(fApprox)
% loglog(Ns,PowerBound+delta,'k--','linewidth',1.5)
% 
% xlabel('$N$','interpreter','latex')
% ylabel('$\|\tilde{f}_N\|$','interpreter','latex')
% set(gca,'fontsize',20)
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% legend('After 10 s','After 100 s','After 1000 s','After 2000 s')
% %% Post Processing
% % Theta_hat_E=simresult.Theta_E.Data;
% Theta_hat_3m=simresult.Theta_3m.Data;
% Theta_hat_5m=simresult.Theta_5m.Data;
% tsim=simresult.tout;
% 
