
delta = 0;
triggered = 1;
ext = 0;
A=-1e4; B=1;
Ar=-1e4; Br=1;
ga = 0; 
gb = 0;
gf = 1e7;
fMag = 2;
SNU = 2; %Signal to Uncertainty
OmegaLim = 1.75;
Scl = 0.7;


x0=OmegaLim-1.5; xr0=OmegaLim-0.05;

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

runTimes = 10;
Newt = 1;
model = 'NewtTrigger';
open(model)
trackError0 = abs(x0-xr0);
maxN = 90;
C = 1;


ssErrorNB = [];
Epsilons = [];

cs = getActiveConfigSet(model);
csNew = copy(cs);
set_param(csNew,"StopTime",num2str(runTimes))
StartMeasureTime = runTimes-10;
startTime = 0;
h = 1e-5;
N0 = 0;
lambda0 = [];
triggerErrors = [];
triggerTimes = [];
trackingError = [];
totTime = [];
tol = C*trackError0;
while trackError0 > 1e-5 && N0 < maxN && startTime < runTimes
set_param(csNew,"StartTime",num2str(startTime))


[gCenters,offlineError,powBound,~,~,~,~,gV] = greedyCenters(tol,h,OmegaLim,ext,maxN,type,par,scale,interps,targetType,targetNu,targetScale,targetCoefs,Dim,fMag,Scl);
norm_fs = 4.42;
if size(gCenters,1) > N0
    C = 0.5*C;
    N = size(gCenters,1)
    lambda0 = [lambda0;zeros(N-N0,1)];
    if N > 1
        PowerBound = powBound(N+1)
    else
        PowerBound = powBound(N)
    end
    tol = C*PowerBound;
    N0 = N;
    gK = kermat(gCenters,gCenters,type,par,scale);
    
    [L,D,Perm] = ldl(gK);
    
    
    L2 = chol(gK);
    l = 1./diag(L2);
    L2 = L2*diag(l);
    
    gV = gV(1:N,1:N);
    z = diag(gK);
    P = sqrt(gK(1,1));
    
    
    V = zeros(N,N);
    
    V(:,1) = gK(:,1)/P;
    if size(V,1) > 1
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
    end
end


epsilon = PowerBound+delta;

Epsilons = [Epsilons;epsilon];

%% RKHS Preparation
centers=gCenters;
K = kermat(centers,centers,type,par,scale);

%% Run Simulink
% 
simresult=sim(model,csNew);


xNB = simresult.xNB.Data;
time = simresult.xNB.time;


xr = simresult.xr.Data;

fNB = simresult.fNB.Data;
fhatNB = -simresult.fhatNB.Data;

nonCo = simresult.nonCo.Data;

controller = simresult.inputNB.Data;
lambda = simresult.lambdaHat.Data;

r = simresult.reference.Data;
x0 = xNB(end);
xr0 = xr(end);

lambda0 = lambda(end,:)';
trackError0 = abs(simresult.trackingErrorNB.Data(end))
trackingError = [trackingError;abs(simresult.trackingErrorNB.Data)];
startTime = time(end)
totTime = [totTime;time];
triggerTimes = [triggerTimes;startTime];
triggerErrors = [triggerErrors; trackError0];
end
triggered = 0;
set_param(csNew,"StartTime",num2str(startTime))
% 
% 


%% Run Simulink
% 
simresult=sim(model,csNew);


xNB = simresult.xNB.Data;
time = simresult.xNB.time;


xr = simresult.xr.Data;

fNB = simresult.fNB.Data;
fhatNB = -simresult.fhatNB.Data;

nonCo = simresult.nonCo.Data;

controller = simresult.inputNB.Data;
lambda = simresult.lambdaHat.Data;

r = simresult.reference.Data;

trackingError = [trackingError;abs(simresult.trackingErrorNB.Data)];
totTime = [totTime;time];



%%
fSS = figure(2);

hold on



semilogy(triggerTimes(1:end-1),triggerErrors(1:end-1),'o','linewidth',1.5);

semilogy(totTime,trackingError,'linewidth',1.5);
tm = [0;triggerTimes];
line([tm(end-1),runTimes],[Epsilons(end-1), Epsilons(end-1)],'linestyle','--','linewidth',2,'color','k')
line([tm(end-1),runTimes],[Epsilons(end), Epsilons(end)],'linestyle','-.','linewidth',1.5,'color','k')
for ll = 2:length(triggerTimes)
    line([tm(ll-1),tm(ll)],[Epsilons(ll-1), Epsilons(ll-1)],'linestyle','-.','linewidth',1.5,'color','k')
    line([tm(ll),tm(ll)],[Epsilons(ll-1), Epsilons(ll)],'linestyle','-.','linewidth',1.5,'color','k')
end



% loglog(delta*ones(length(trackErrors),1),'k-.','linewidth',1.5)
xlabel('log(time (s))','interpreter','latex')
ylabel('log($\frac{\|\tilde{x}_N\|}{\|B^TP\|})$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'YScale','log')
set(gca,'XScale','log')
legend('trigger events at $t_m$','tracking error','ultimate error bound','$\sup_{\xi_N \in \Omega N(m)} P_{H_{N}}$','interpreter','latex')
% error = simresult.trackingErrorNB.Data;
% 
% 
% figure(); plot(time,error,'linewidth',1.5); hold on; plot(time,xNB,'linewidth',1.5);
%  plot(time,xr,'linewidth',1.5);
% title(strcat('\gamma = ',num2str(gf)))
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
% 
% 
% 
