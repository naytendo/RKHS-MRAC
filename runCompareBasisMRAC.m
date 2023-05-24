


delta = 0;

ext = 0;
A=-1; B=1;
Ar=-1; Br=1;
ga = 0; 
gb = 0;
gf = 2e-2;
fMag = 2;
SNU = 2; %Signal to Uncertainty
OmegaLim = 1.75;
Scl = 0.7;


x0=OmegaLim-0.05; xr0=OmegaLim-0.05;
alpha0=0; beta0=0;sigma=0.1;
kx0=0; kr0=0;

% type = 'exp';
type = 'g';
Dim = 1;
scale=0.1;

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
ySamps = fMag*tanh(Scl*interps);

targetCoefs = intmat\ySamps;


tol = 1e-13;
maxN = 90;
h = 1e-4;
[gCenters,offlineError,powBound,norm_fs,condNum,condInt,condNew2,gV] = greedyCenters(tol,h,OmegaLim,ext,maxN,type,par,scale,interps,targetType,targetNu,targetScale,targetCoefs,Dim,fMag,Scl);


N = length(gCenters);
gK = kermat(gCenters,gCenters,type,par,scale);

[L,D,Perm] = ldl(gK);


L2 = chol(gK);
l = 1./diag(L2);
L2 = L2*diag(l);

gV = gV(1:N,1:N);
z = diag(gK);
v_x = zeros(N,1);
P = sqrt(gK(1,1));


V = zeros(N,N);

V(:,1) = gK(:,1)/P;
w = V(:,1).^2;
P = sqrt(z(2)-w(2));

for jj = 1:N-1
    e = zeros(N,1);
    e(jj+1) = 1;
    Vm = V(:,1:jj);
    u = gK*e-Vm*Vm'*e;
    V(:,jj+1) = u/P; % Newton Basis evaluated at the centers
    w = V(:,jj+1).^2 + w; % Updating sum of Newton basis squared at the centers
    if jj ~=N-1
        P = sqrt(z(jj+2)-w(jj+2));
    end
end
figure(1)
loglog(condNum,'-s','linewidth',2)
hold on
loglog(condInt,'-o','linewidth',2)

ylabel('Condition Number of Interpolation Matrix')
xlabel('$N$','interpreter','latex','fontsize',20)




runTimes = 50;
Newt = 1;
model = 'newtVsStdMRAC';
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
StartMeasureTime = runTimes(tt)-10;



Ns = floor(logspace(log10(4),log10(length(gCenters)),10));

Ns = Ns(end-4:end);


for nn = Ns
N = nn

epsilon = powBound(nn)*norm_fs+delta;

%% RKHS Preparation
centers=gCenters(1:nn);
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

fSS = figure(2);

hold on
loglog(ssError(Ns),'-s','linewidth',1.5);
loglog(ssErrorNB(Ns),'-o','linewidth',1.5);


% loglog(Ns,PowerBound,'k-.','linewidth',1.5);


% fApprox = figure(12);
% hold on
% loglog(Ns,approxError,'-o','linewidth',1.5)
% error = simresult.trackingErrorNB.Data;
% xNB = simresult.xNB.Data;
% 
% xr = simresult.xr.Data;
% time = simresult.ssError.time;
% 
% fNB = simresult.fNB.Data;
% fhatNB = -simresult.fhatNB.Data;
% 
% nonCo = simresult.nonCo.Data;
% 
% controller = simresult.inputNB.Data;
% 
% r = simresult.reference.Data;
% 
% % figure(); plot(time,error,'linewidth',1.5); hold on; plot(time,xNB,'linewidth',1.5);
% %  plot(time,xr,'linewidth',1.5);
% % title(strcat('\gamma = ',num2str(gf)))
% 
% 
% plot(time,(powBound(end)*norm_fs+delta)*ones(length(time),1),'k--','linewidth',1.5)
% plot(time,-(powBound(end)*norm_fs+delta)*ones(length(time),1),'k--','linewidth',1.5)
% legend('$\tilde{x}_N$','$x_N$','$x_r$','$\epsilon$','interpreter','latex')
% 
% 
% figure(); plot(time,fNB,'linewidth',1.5); hold on;
% plot(time,-nonCo,'linewidth',1.5);plot(time,controller,'linewidth',1.5);plot(time,r,'linewidth',1.5);
% title(strcat('\gamma = ',num2str(gf)))
% 
% legend('$f$','$\hat{f}_N$','$u_N$','$r$','$\epsilon$','interpreter','latex')
% 
% 
% figure(); plot(time,fNB,'linewidth',1.5); hold on; plot(time,fhatNB,'linewidth',1.5)
% 
% legend('$f$','$\hat{f}_N$','interpreter','latex')
% title(strcat('\gamma = ',num2str(gf)))
% % loglog(Ns,PowerBound,'k-.','linewidth',1.5);

%%
figure(fSS)
loglog(powBound(Ns)*norm_fs+delta,'k--','linewidth',1.5);

loglog(delta*ones(length(Ns),1),'k-.','linewidth',1.5)
xlabel('$N$','interpreter','latex')
ylabel('$\frac{\|\tilde{x}_N\|}{\|B^TP\|}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')
legend('$\gamma$ = 1e-2','$\gamma$ = 1e-3','$\gamma$ = 1e-4','$\gamma$ = 1e-5','$\epsilon$','interpreter','latex')
