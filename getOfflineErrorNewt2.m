function [offlineErrorNew,PowerBound,condNum,norm_fs,fineCenters] = getOfflineErrorNewt2(Ns,type,RBFpar,scale,OmegaLim,targetNu,RBFtargetType,targetScale,targetCoefs,interps,Dim,Newt,ext)

% Dim = 1;
% % hype=1/(2*2^2); 
% hyperMtrn=2; sig=2; 
% 
% nu = 1/2*(10+Dim/2);
% %% Target Function generation
% s = 8;
% targetNu = 1/2*(s+Dim/2);
% Ip = 8;
% targetHyper = 1;
% Ns = [2+1,4+1,8+1,16+1,32+1,64+1,128+1];
% Ns = 2:128;
Nmesh = 1024+1;
xMesh = linspace(-OmegaLim,OmegaLim,Nmesh)';
% KMesh = kermat(xMesh,interps,RBFtargetType,targetNu-Dim/2,targetScale);
% yMesh = KMesh*targetCoefs;
yMesh = tanh(0.1*xMesh);

norm_fs = zeros(length(Ns),1);
norm_fsNew = zeros(length(Ns),1);

% Ns = [4+1,8+1,16+1,32+1,64+1,128+1];

condNum = zeros(length(Ns),1);
condNumNew = zeros(length(Ns),1);
offlineError = zeros(length(Ns),1);
offlineErrorNew = zeros(length(Ns),1);
PowerBound = zeros(length(Ns),1);


for nn = 1:length(Ns)
    centers = linspace(-OmegaLim-ext,OmegaLim+ext,Ns(nn))';
    N = length(centers);

    Kc = kermat(centers,centers,type,RBFpar,scale);
    KcMesh = kermat(centers,xMesh,type,RBFpar,scale);
%     U=Kc\KcMesh;
%     coefs = zeros(max(Ns),1);
%     for ii = 1:length(indexes)   
%         ur(ii,:)=U(ii,:);
%         coefs(ii) = yMesh(indexes(ii))-coefs'*ur(:,indexes(ii));
%     end
%     condNumNewt(nn)=cond(ur(1:Ns(nn),indexes));

    

%     KEval = kermat(centers,interps,RBFtargetType,targetNu-Dim/2,targetScale);
%     yEval = KEval*targetCoefs;

    yEval = tanh(0.1*centers);
    
    condNum(nn) = cond(Kc);
    theta =Kc\yEval;


    V = zeros(N,N);
    z = diag(Kc);

    P = sqrt(real(kermat(centers(1),centers(1),type,RBFpar,scale)));

    V(:,1) = real(kermat(centers(1),centers,type,RBFpar,scale))/P;
    w = V(:,1).^2;
    P = sqrt(z(2)-w(2));

    for jj = 1:N-1
        e = zeros(N,1);
        e(jj+1) = 1;
        Vm = V(:,1:jj);
        u = real(kermat(centers,centers(jj+1),type,RBFpar,scale))-Vm*Vm'*e;
        V(:,jj+1) = u/P;
        w = V(:,jj+1).^2 + w;
        if jj < N-1
            P = sqrt(z(jj+2)-w(jj+2));
        end
        
    end

    coefsNew = zeros(N,1);
    coefsNew(1) = yEval(1)/V(1,1);
    for kk = 2:N
        coefsNew(kk) = yEval(kk)/V(kk,kk)-V(kk,1:kk-1)*coefsNew(1:kk-1)/V(kk,kk);
    end

%     coefsNew = V\yEval;
    condNumNew(nn) = max(diag(V))/min(diag(V));



%     norm_fs(index) = sqrt(coefs(1:Ns(nn))'*ur(1:Ns(nn),indexes)*coefs(1:Ns(nn)));
%    
    norm_fs(nn) = sqrt(theta'*Kc*theta);
    norm_fsNew(nn) = sqrt(coefsNew'*coefsNew);
    y = KcMesh'*theta;
    
    yNew = zeros(length(xMesh),1);
    for mm = 1:length(xMesh)
        nbNew = newtonBasisAtX(xMesh(mm),centers,type,RBFpar,scale);
        yNew(mm) = coefsNew'*nbNew;
    end

%     err=yEst'-yMesh;
    
    err = y - yMesh;
    errNew = yNew - yMesh;

%     offlineError(nn)=max(abs(err));
    offlineError(nn)=max(abs(err));
    offlineErrorNew(nn)=max(abs(errNew));

    
    
end

fineCenters = centers;
figure(8)
hold on;
% plot(Ns,norm_fs,'-o','linewidth',2)
plot(Ns,norm_fs,'-s','linewidth',2)
plot(Ns,norm_fsNew,'-x','linewidth',2)
xlabel('$N$','interpreter','latex')
ylabel('$\| \Pi_N f \|_{\mathcal{H}}$','interpreter','latex')
plot(Ns(condNum < 10^5),max(norm_fs(condNum < 10^5))*ones(length(Ns(condNum < 10^5)),1),'k--','linewidth',1.5)
legend('Standard Basis','Newton Basis','$\approx \|f\|_{\mathcal{H}}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

% NMax = max(NMax);

norm_fs = max(norm_fs(condNum < 10^5));
% Ns = [4+1,8+1,16+1,32+1,64+1,128+1];

figure(1)
hold on
plot(xMesh,yMesh,'linewidth',4,'color',[0 0.4470 0.7410])
plot(xMesh,yNew,'--','linewidth',2,'color',[0.9290 0.6940 0.1250])

% title('verify good approximation of matched uncertainty')

legend('Actual','Approximation','location','southeast')
xlabel('\Omega')
ylabel('f(x)')
set(gca,'fontsize',20)

% xEval = linspace(-OmegaLim,OmegaLim,1028+1)';


% 
figure(3)
hold on
loglog(Ns,condNum,'-o','linewidth',1.5)
if Newt
    loglog(Ns,condNumNew,'-d','linewidth',1.5)
    legend('Standard Basis','Newton Basis')
end
xlabel('$N$','interpreter','latex')
ylabel('Cond($A_\psi(\Xi_N)$)','interpreter','latex')

set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

for nn = 1:length(Ns)
    centers = linspace(-OmegaLim-ext,OmegaLim+ext,Ns(nn))';
    Kc = kermat(centers,centers,type,RBFpar,scale);
    Pow = zeros(length(xMesh),1);
    for mm = 1:length(xMesh)
        kv = kermat(centers,xMesh(mm),type,RBFpar,scale);
        Pow(mm) = powerFunction(type,xMesh(mm),kv,Kc,scale,RBFpar);
    end
    PowerBound(nn) = max(Pow)*norm_fs;

end

figure(2)
hold on

loglog(Ns,offlineError,'-o','linewidth',1);
loglog(Ns,offlineErrorNew,'-d','linewidth',1);
% colorApprox = get(lineApprox,'Color');
loglog(Ns,PowerBound,'k-.','linewidth',1.5);
ylabel('$\|f\|_{\mathcal{H}} - \|\Pi_N f\|_{\mathcal{H}}$','interpreter','latex')
legend('Standard Basis','Newton Basis','$\sup_{x \in \Omega}P_{\mathcal{H}_N}(x)\|f\|_{\mathcal{H}}$','interpreter','latex')
xlabel('$N$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')





end

