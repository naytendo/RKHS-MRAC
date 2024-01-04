Ri = 2*sqrt(2)-0.3;
Ro = 2*sqrt(2)+0.3;
delta = 0;
triggered = 1;
del = 0;

A=[0 1;-1 -1]; B=[0;1];
Ar=A; Br=B;
x0=[2,3]; xr0=[2,2];
rmax = norm(xr0);
gf = 1e5;

%% verify controllability
if  rank(ctrb(A,B)) == size(A,1)
    fprintf('System controllable\n')
    pDes = [-0.5000 + 0.8660i;-0.5000 - 0.8660i];
    kx = place(A,B,pDes);
else
    fprintf('System uncontrollable\n')
end
ga = 0; 
gb = 0;

Q = eye(2);

P = lyap(A,Q);

mu = 0.01;





% type = 'exp';
type = 'ms';
scale=0.3;
par = 6;

%% Target Function generation



runTimes = 2000;
Newt = 1;
model = 'NewtTriggerVanDerPol';
open(model)
trackError0 = abs(x0-xr0);
maxN = 70;

ssErrorNB = [];


cs = getActiveConfigSet(model);
csNew = copy(cs);
set_param(csNew,"StopTime",num2str(runTimes))
StartMeasureTime = runTimes-10;
startTime = 0;
h = 1e-1;
N0 = 0;
lambda0 = [];
triggerErrors = [];
triggerTimes = [];
trackingError = [];
Epsilons = [];
totTime = [];
rTot = [];
xTot = [];
xrTot = [];
input = [];
vn = [];
uSL = [];
nonCo = [];
fhat = [];
fNon = [];
Ns = [];
condK = [];
condN =[];
C = 0.8;

Omega = omegaCover2(Ri,Ro,h);
%%
[norm_fs,condT,Ntest] = getVDP_HNorm(type,par,scale,mu,Ri,Ro);

norm_f = max(norm_fs(condT < 1e12))
format short g;
condT
tol = C*norm(x0-xr0);
R = norm_f+0.1;
triggerTime = 0;

%%
while norm(trackError0) > 1e-7 && N0 < maxN && triggerTime < runTimes
    set_param(csNew,"StartTime",num2str(triggerTime))
    
    
    [gCenters,offlineError,powBound,condNum,condNumNewt,gV] = greedySet2(tol,Omega,maxN,type,par,scale,norm_f);
    
    N = size(gCenters,1)
    Ns = [Ns;N];
    
    if N > N0
        
        N = size(gCenters,1);
        lambda0 = [lambda0;zeros(N-N0,1)];
        
        if N > 1 && N ~= maxN
            PowerBound = powBound(N+1);
        else
            PowerBound = powBound(N);
        end
        N0 = N;
        gK = real(kermat(gCenters,gCenters,type,par,scale));
        
        
        
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
                uSL = gK*e-Vm*Vm'*e;
                V(:,jj+1) = uSL/P; % Newton Basis evaluated at the centers
                w = V(:,jj+1).^2 + w; % Updating sum of Newton basis squared at the centers
                if jj ~=N-1
                    P = sqrt(z(jj+2)-w(jj+2));
                end
            end
        end
    end
    tol = C*PowerBound*R
    condN = [condN;cond(V)];
    epsilon = PowerBound*R+delta
    
    Epsilons = [Epsilons;epsilon];
    
    %% RKHS Preparation
    centers=gCenters;
    K = real(kermat(centers,centers,type,par,scale));
    condK = [condK;cond(K)];
    %% Run Simulink
    % 
    simresult=sim(model,csNew);
    
    %%
    x = simresult.x.Data;
    time = simresult.x.time;
    triggerTime = time(end)
    
    xr = simresult.xr.Data;
    
    fNon = [fNon;simresult.fNon.Data];
    fhat = [fhat;-simresult.fhat.Data];
    
    nonCo = [nonCo;simresult.nonCo.Data];
    vn = [vn;simresult.vn.Data];
    uSL = [uSL;simresult.uSL.Data];
    
    input = [input;simresult.u.Data];
    lambda = simresult.lambdaHat.Data;
    
    r = simresult.reference.Data;
    x0 = x(end,:);
    xr0 = xr(end,:);
    rTot = [rTot;r];
    xrTot = [xrTot;xr];
    xTot = [xTot;x];
    lambda0 = lambda(end,:)';
    trackError0 = norm(simresult.trackingError.Data(end,:))
    trackingError = [trackingError;abs(simresult.trackingError.Data)];
    totTime = [totTime;time];
    if triggerTime < runTimes && N0 < maxN && trackError0 > 1e-7
        triggerErrors = [triggerErrors;trackError0];
        triggerTimes = [triggerTimes;triggerTime];
        
    end
    if N0 == maxN
        triggered = 0;
    end
end
triggered = 0;
if triggerTime < runTimes
    set_param(csNew,"StartTime",num2str(triggerTime))
    % 
    % 
    
    
    %% Run Simulink
    % 
    simresult=sim(model,csNew);
    
    %%
    x = simresult.x.Data;
    time = simresult.x.time;
    
    
    xr = simresult.xr.Data;
    
    fNon = [fNon;simresult.fNon.Data];
    fhat = [fhat;-simresult.fhat.Data];
    
    nonCo = [nonCo;simresult.nonCo.Data];
    vn = [vn;simresult.vn.Data];
    uSL = [uSL;simresult.vn.Data];
    
    input = [input;simresult.u.Data];
    lambda = simresult.lambdaHat.Data;
    
    r = simresult.reference.Data;
    
    trackingError = [trackingError;abs(simresult.trackingError.Data)];
    totTime = [totTime;time];
    xTot = [xTot;x];
    xrTot = [xrTot;xr];
end


%%
fSS = figure();



subplot(3,1,1:2)
hold on
% trackingErrorNorm = sqrt(sum(trackingError.^2,2));
semilogy(totTime,trackingErrorNorm,'linewidth',1.5);




if ~isempty(triggerTimes)
    tm = [totTime(2);triggerTimes];
    line([tm(1),runTimes],[Epsilons(end), Epsilons(end)],'linestyle','--','linewidth',2.5,'color',[0.9290 0.6940 0.1250])
    l = semilogy(triggerTimes(1:end),triggerErrors(1:end),'o','markerSize',5);
    l.MarkerFaceColor = l.Color;
    line([tm(end),runTimes],[Epsilons(end), Epsilons(end)],'linestyle','-.','linewidth',1.5,'color','k')
    for ll = 1:length(triggerTimes)
        line([tm(ll),tm(ll+1)],[Epsilons(ll), Epsilons(ll)],'linestyle','-.','linewidth',1.5,'color','k')
        line([tm(ll+1),tm(ll+1)],[Epsilons(ll), Epsilons(ll+1)],'linestyle','-.','linewidth',1.5,'color','k')
    end
end


set(gca,'Yscale','log')
set(gca,'fontsize',10)
ylabel('$\frac{\|\tilde{x}_N\|}{\|B^TP\|}$','interpreter','latex','fontsize',20)

legend('tracking error','ultimate error bound','trigger events at $t_m$','$\sup_{\xi \in \Omega_{m}} P_{H_{N}}(\xi)$','interpreter','latex','fontsize',20)
grid on
subplot(3,1,3)
plot([tm;runTimes],[Ns;Ns(end)],'linewidth',1.5)
xlabel('time (s)','interpreter','latex','fontsize',20)
ylabel('$N$','interpreter','latex','fontsize',20)
ylim([10 500])
set(gca,'Yscale','log')
grid on
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
%%
% figure(); plot(totTime,fNon,'linewidth',2); hold on;
% plot(totTime,-nonCo,'linewidth',1); plot(totTime,vn,'linewidth',2);
% plot(totTime,input,'linewidth',1);
% title(strcat('\gamma = ',num2str(gf)))
% 
% legend('$f$','$-\hat{f}_N$','$v_N$','$u_{SL}$','$u_N$','interpreter','latex')
% % 

%%
% 
% figure(); plot(time,fNon,'linewidth',1.5); hold on; plot(time,fhat,'linewidth',1.5)
% 
% legend('$f$','$\hat{f}_N$','interpreter','latex')
% title(strcat('\gamma = ',num2str(gf)))
% loglog(Ns,PowerBound,'k-.','linewidth',1.5);
% 
% 
% 
%%
figure();
% plot(Omega(:,1),Omega(:,2),'k.','markersize',1)
hold on
plot(gCenters(:,1),gCenters(:,2),'ko')
plot(xTot(:,1),xTot(:,2),'b-.','linewidth',2)
grid on
theta = linspace(0,360,100);
plot(rmax*sind(theta),rmax*cosd(theta),'color',[0.4940 0.1840 0.5560],'linewidth',5)
xlabel('$x_1$','interpreter','latex','fontsize',20)
ylabel('$x_2$','interpreter','latex','fontsize',20)


index = find(totTime == triggerTimes(end-1));
plot(xTot(index:end,1),xTot(index:end,2),'color',[0.9290 0.6940 0.1250],'linewidth',1.5)
legend('$\Xi_N$','$x$','$x_r$','$x_{ss}$','fontsize',20,'interpreter','latex')