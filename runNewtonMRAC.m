

delta = 0;
ext = 0.1;
A=-4; B=0.8;
Ar=-5; Br=1;
ga = 10; 
gb = 100; 
gf = 10;
fMag = 1;
SNU = 10; %Signal to Uncertainty
OmegaLim = 7;


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
tol = 0;
maxN = 128;
h = 1e-3;
[gCenters,error,powBound,norm_fs,condNum,condNew] = greedyCenters(tol,h,OmegaLim,ext,maxN,type,par,scale,interps,targetType,targetNu,targetScale,targetCoefs,Dim);
% figure()
% plot([-OmegaLim-ext, OmegaLim+ext],[0 0],'k','linewidth',2);
% hold on
% showCenter = 10;
% plot(gCenters(1:showCenter),0*gCenters(1:showCenter),'ro','linewidth',2);
% for cc = 1:showCenter
%     text(gCenters(cc)-0.1,-0.005,num2str(cc),'fontsize',20);
% end
% xlim([-OmegaLim-ext,OmegaLim+ext])
% ylim([-0.04 0.04]);


figure()
loglog(condNum,'-s','linewidth',2)
hold on
loglog(condNew,'-o','linewidth',2)

ylabel('Condition Number of Interpolation Matrix')
xlabel('$N$','interpreter','latex','fontsize',20)


Ns = [4+1,8+1,16+1,32+1,64+1];

runTimes = 50;
Newt = 1;

[offlineError,PowerBound,condNum,norm_fs] = getOfflineErrorNewt2(Ns,type,par,scale,OmegaLim,targetNu,targetType,targetScale,targetCoefs,interps,Dim,Newt,ext);

model = 'NewtonMRAC';
open(model)
for tt = 1:length(runTimes)
ssErrorWendland = zeros(length(gCenters),1);
approxError = zeros(length(gCenters),1);
cs = getActiveConfigSet(model);
csNew = copy(cs);
set_param(csNew,"StopTime",num2str(runTimes(tt)))
StartMeasureTime = runTimes(tt)-10;


for nn = 1:10:64
N = nn;

epsilon = powBound(nn)*norm_fs+delta;

%% RKHS Preparation
centers=gCenters(1:nn);
K = kermat(centers,centers,type,par,scale);
% condNum(nn)=cond(K);
% L = chol(K,'lower');
% condNumKL(nn) = max(diag(L))/min(diag(L));





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
lineSS = loglog(1:10:64,ssErrorWendland(1:10:64),'-d','linewidth',1.5);
colorSS = get(lineSS,'color');

% loglog(Ns,PowerBound,'k-.','linewidth',1.5);


% fApprox = figure(12);
% hold on
% loglog(Ns,approxError,'-o','linewidth',1.5)


% loglog(Ns,PowerBound,'k-.','linewidth',1.5);

end
%%
figure(fSS)
loglog(1:10:64,powBound(1:10:64)*norm_fs+delta,'k--','linewidth',1.5);

% loglog(Ns,delta*ones(length(Ns),1),'k-.','linewidth',1.5)
xlabel('$N$','interpreter','latex')
ylabel('$\frac{\|\tilde{x}_N\|}{\|B^TP\|}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('Standard Basis','Newton Basis','$\epsilon$','interpreter','latex')
% legend('\gamma_b = 10','\gamma_b= 100','\gamma_b = 500','\gamma_b = 1000','\epsilon')

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
