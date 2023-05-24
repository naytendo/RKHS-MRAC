function [offlineErrorRec,PowerBound,condNum,norm_fs] = getOfflineErrorNewt(Ns,type,RBFpar,scale,OmegaLim,targetNu,RBFtargetType,targetScale,targetCoefs,interps,Dim,Newt,ext)

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
KMesh = kermat(xMesh,interps,RBFtargetType,targetNu-Dim/2,targetScale);
yMesh = KMesh*targetCoefs;


norm_fsStand = zeros(length(Ns),1);
norm_fsRec = zeros(length(Ns),1);
norm_fsDir = zeros(length(Ns),1);

% Ns = [4+1,8+1,16+1,32+1,64+1,128+1];

condNum = zeros(length(Ns),1);
condNumRec = zeros(length(Ns),1);
condNumDir = zeros(length(Ns),1);
offlineErrorStand = zeros(length(Ns),1);
offlineErrorDir = zeros(length(Ns),1);
offlineErrorRec = zeros(length(Ns),1);
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

    [R,P] = pivchol(Kc);
 

    L = P*R';

    L2 = chol(Kc,'lower');

    condNumDir(nn) = cond(L);

    KEval = kermat(centers,interps,RBFtargetType,targetNu-Dim/2,targetScale);
    yEval = KEval*targetCoefs;
    
    condNum(nn) = cond(Kc);
    theta =Kc\yEval;
    coefsDir = L\yEval;

    coefsRec = zeros(N,1);
    coefsRec(1) = yEval(1)/L(1,1);
    for kk = 2:N
        coefsRec(kk) = yEval(kk)/L(kk,kk)-L(kk,1:kk-1)*coefsRec(1:kk-1)/L(kk,kk);
    end

    condNumRec(nn) = max(diag(L))/min(diag(L));

%     norm_fs(index) = sqrt(coefs(1:Ns(nn))'*ur(1:Ns(nn),indexes)*coefs(1:Ns(nn)));
%    
    norm_fsStand(nn) = sqrt(theta'*Kc*theta);
    norm_fsDir(nn) = sqrt(coefsDir'*(L\Kc/L')*coefsDir);
    norm_fsRec(nn) = sqrt(coefsRec'*coefsRec);


    
%     yEst=coefs'*ur;
    yDir = (L\KcMesh)'*coefsDir;
    yStand = KcMesh'*theta;

    yRec = zeros(length(xMesh),1);

    for mm = 1:length(xMesh)
        knlV = kermat(centers,xMesh(mm),type,RBFpar,scale);
        nb = zeros(N,1);
        nb(1) = 1/L(1,1)*knlV(1);
        for kk = 2:N
            nb(kk) = knlV(kk)/L(kk,kk)-L(kk,1:kk-1)*nb(1:kk-1)/L(kk,kk);
        end
        yRec(mm) = coefsRec'*nb;
    end

%     err=yEst'-yMesh;
    
    errStand = yStand - yMesh;
    errDir = yDir - yMesh;
    errRec = yRec - yMesh;

%     offlineError(nn)=max(abs(err));
    offlineErrorStand(nn)=max(abs(errStand));
    offlineErrorDir(nn)=max(abs(errDir));
    offlineErrorRec(nn)=max(abs(errRec));

    
    
end


figure(8)
hold on;
% plot(Ns,norm_fs,'-o','linewidth',2)
plot(Ns,norm_fsStand,'-s','linewidth',2)
plot(Ns,norm_fsDir,'-d','linewidth',2)
plot(Ns,norm_fsRec,'-o','linewidth',2)
xlabel('$N$','interpreter','latex')
ylabel('$\| \Pi_N f \|_{\mathcal{H}}$','interpreter','latex')
plot(Ns(condNum < 10^5),max(norm_fsDir(condNum < 10^5))*ones(length(Ns(condNum < 10^5)),1),'k--','linewidth',1.5)
legend('Standard Basis','Newton Basis','Newton Basis Recursive','$\approx \|f\|_{\mathcal{H}}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

% NMax = max(NMax);

norm_fs = max(norm_fsStand(condNum < 10^5));
% Ns = [4+1,8+1,16+1,32+1,64+1,128+1];

figure(1)
hold on
plot(xMesh,yMesh,'linewidth',4,'color',[0 0.4470 0.7410])
plot(xMesh,yRec,'--','linewidth',2,'color',[0.9290 0.6940 0.1250])

% title('verify good approximation of matched uncertainty')

legend('Actual','Approximation','location','southeast')
xlabel('\Omega')
ylabel('f(x)')
set(gca,'fontsize',20)

% xEval = linspace(-OmegaLim,OmegaLim,1028+1)';


% 
figure(3)
hold on
loglog(Ns,condNum,'-s','linewidth',1.5)
if Newt
    loglog(Ns,condNumDir,'-d','linewidth',1.5)
    loglog(Ns,condNumRec,'-o','linewidth',1.5)
    legend('Standard Basis','Newton Basis Direct','Newton Basis Recursion')
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

loglog(Ns,offlineErrorStand,'-o','linewidth',1);
loglog(Ns,offlineErrorDir,'-s','linewidth',1);
loglog(Ns,offlineErrorRec,'-d','linewidth',1);
% colorApprox = get(lineApprox,'Color');
loglog(Ns,PowerBound,'k-.','linewidth',1.5);
ylabel('$\|f\|_{\mathcal{H}} - \|\Pi_N f\|_{\mathcal{H}}$','interpreter','latex')
legend('Standard Basis','Newton Basis Direct','Newton Basis Recursion','$\sup_{x \in \Omega}P_{\mathcal{H}_N}(x)\|f\|_{\mathcal{H}}$','interpreter','latex')
xlabel('$N$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')





end

