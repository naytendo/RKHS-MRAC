function [offlineError,PowerBound,condNum,norm_fs] = getOfflineError(Ns,type,RBFpar,scale,OmegaLim,targetNu,RBFtargetType,targetScale,targetCoefs,interps,Dim,ext)
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


norm_fs = zeros(length(Ns),1);


% Ns = [4+1,8+1,16+1,32+1,64+1,128+1];

condNum = zeros(length(Ns),1);
offlineError = zeros(length(Ns),1);
PowerBound = zeros(length(Ns),1);


for nn = 1:length(Ns)
    centers = linspace(-OmegaLim-ext,OmegaLim+ext,Ns(nn))';

    Kc = kermat(centers,centers,type,RBFpar,scale);
    KcMesh = kermat(centers,xMesh,type,RBFpar,scale);
%     U=Kc\KcMesh;
%     coefs = zeros(max(Ns),1);
%     for ii = 1:length(indexes)   
%         ur(ii,:)=U(ii,:);
%         coefs(ii) = yMesh(indexes(ii))-coefs'*ur(:,indexes(ii));
%     end
%     condNumNewt(nn)=cond(ur(1:Ns(nn),indexes));


    KEval = kermat(centers,interps,RBFtargetType,targetNu-Dim/2,targetScale);
    yEval = KEval*targetCoefs;
    
    condNum(nn) = cond(Kc);
    theta =Kc\yEval;


%     norm_fs(index) = sqrt(coefs(1:Ns(nn))'*ur(1:Ns(nn),indexes)*coefs(1:Ns(nn)));
%    
    norm_fs(nn) = sqrt(theta'*Kc*theta);



    
%     yEst=coefs'*ur;
    yEst = KcMesh'*theta;
    err = yEst - yMesh;
    offlineError(nn)=max(abs(err));


    
    
end


figure(8)
hold on;
% plot(Ns,norm_fs,'-o','linewidth',2)
plot(Ns,norm_fs,'-s','linewidth',2)

xlabel('$N$','interpreter','latex')
ylabel('$\| \Pi_N f \|_{\mathcal{H}}$','interpreter','latex')
plot(Ns(condNum < 10^5),max(norm_fs(condNum < 10^5))*ones(length(Ns(condNum < 10^5)),1),'k--','linewidth',1.5)
% legend('Standard Basis','Newton Basis','Newton Basis Recursive','$\approx \|f\|_{\mathcal{H}}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

% NMax = max(NMax);

norm_fs = max(norm_fs(condNum < 10^5));
% Ns = [4+1,8+1,16+1,32+1,64+1,128+1];

figure(1)
hold on
plot(xMesh,yMesh,'linewidth',4,'color',[0 0.4470 0.7410])
plot(xMesh,yEst,'--','linewidth',2,'color',[0.9290 0.6940 0.1250])

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

% colorApprox = get(lineApprox,'Color');
loglog(Ns,PowerBound,'k-.','linewidth',1.5);
ylabel('$\|f\|_{\mathcal{H}} - \|\Pi_N f\|_{\mathcal{H}}$','interpreter','latex')
% legend('Standard Basis','Newton Basis Direct','Newton Basis Recursion','$\sup_{x \in \Omega}P_{\mathcal{H}_N}(x)\|f\|_{\mathcal{H}}$','interpreter','latex')
xlabel('$N$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

end
