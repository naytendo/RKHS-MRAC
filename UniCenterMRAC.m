


delta = 0.005;

ext = 0;
A=-1; B=1;
Ar=-1; Br=1;
ga = 0; 
gb = 0; 
gf = 1;
fMag = 1;
SNU = 10; %Signal to Uncertainty
OmegaLim = 5;


x0=OmegaLim+2; xr0=OmegaLim+2;
alpha0=0; beta0=0;sigma=0.1;
kx0=0; kr0=0;

% type = 'exp';
type = 'g';
Dim = 1;
scale=0.2;

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


Ns = floor(logspace(log10(4),log10(length(gCenters)),10));
[offlineError,powBound,condNum,norm_fs,uniCenters] = getOfflineErrorNewt2(Ns,type,par,scale,OmegaLim,targetNu,targetType,targetScale,targetCoefs,interps,Dim,Newt,ext);


fineN = length(uniCenters);
uniK = kermat(uniCenters,uniCenters,type,par,scale);
z = diag(uniK);
V = zeros(fineN,fineN);
P = sqrt(uniK(1,1));
V(:,1) = uniK(:,1)/P;
w = V(:,1).^2;
P = sqrt(z(2)-w(2));

for jj = 1:fineN-1
    e = zeros(fineN,1);
    e(jj+1) = 1;
    Vm = V(:,1:jj);
    u = uniK*e-Vm*Vm'*e;
    V(:,jj+1) = u/P; % Newton Basis evaluated at the centers
    w = V(:,jj+1).^2 + w; % Updating sum of Newton basis squared at the centers
    if jj ~=fineN-1
        P = sqrt(z(jj+2)-w(jj+2));
    end
end
runTimes = 30;
Newt = 1;
model = 'CompareBasisMRAC';
open(model)
for tt = 1:length(runTimes)

ssError = zeros(length(gCenters),1);
ssErrorNB = zeros(length(gCenters),1);
approxError = zeros(length(gCenters),1);
condNum = zeros(length(gCenters),1);
condNumNB = zeros(length(gCenters),1);
cs = getActiveConfigSet(model);
csNew = copy(cs);
set_param(csNew,"StopTime",num2str(runTimes(tt)))
StartMeasureTime = runTimes(tt)-1;







for nn = 1:length(Ns)
N = Ns(nn)
centers = linspace(-OmegaLim-ext,OmegaLim+ext,N)';
epsilon = powBound(nn)*norm_fs+delta;

%% RKHS Preparation
K = kermat(centers,centers,type,par,scale);
% L2 = chol(K);
% l = 1./diag(L);
% L2 = L*diag(l);




% 
% 


%% Run Simulink
% 
simresult=sim(model,csNew);
timeSeriesSSError = simresult.ssError;
ssError(nn) = max(timeSeriesSSError.Data);


timeSeriesSSErrorNB = simresult.ssErrorNB;
ssErrorNB(nn) = max(timeSeriesSSErrorNB.Data);

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

fSS = figure(8);

hold on
loglog(Ns,ssError(Ns),'-s','linewidth',1.5);
loglog(Ns,ssErrorNB(Ns),'-o','linewidth',1.5);


% loglog(Ns,PowerBound,'k-.','linewidth',1.5);


% fApprox = figure(12);
% hold on
% loglog(Ns,approxError,'-o','linewidth',1.5)


% loglog(Ns,PowerBound,'k-.','linewidth',1.5);

end
%%
figure(fSS)
loglog(Ns,powBound*norm_fs+delta,'k--','linewidth',1.5);

loglog(Ns,delta*ones(length(Ns),1),'k-.','linewidth',1.5)
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
