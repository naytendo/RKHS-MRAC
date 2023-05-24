function [Xi,error,maxP,norm_fs,condNum,condInt,condNew2,V] = greedySet(tol,Omega,maxN,type,par,scale,interps,targetType,targetNu,targetScale,targetCoefs,Dim,fMag,Scl)




% [Z1,Z2]=meshgrid(-3:h:3,-3:h:3);
% Z=[Z1(:),Z2(:)];


Zmesh = Omega;
% KMesh = kermat(Zmesh,interps,targetType,targetNu-Dim/2,targetScale);
% yMesh = KMesh*targetCoefs;
mu = 3;
yMesh = -(1+mu)-mu*Zmesh(:,1).^2.*Zmesh(:,2);

len=size(Zmesh,1);

% Startmenge [X,Y]
I=floor(len/2);
Xi=Zmesh(I,:);
lambda=zeros(maxN,1);
ur=zeros(maxN,size(Zmesh,1));
error=zeros(maxN,1);
norm_fs=zeros(maxN,1);
maxP=zeros(maxN,1);
intpolnorm=zeros(maxN,1);
newtonnorm=zeros(maxN,1);
condNum=zeros(maxN,1);
condInt=zeros(maxN,1);
condNew2=zeros(maxN,1);
index=[];


for nn=1:maxN
    % Berechnung von A_xz und A_xx
    KXZ=kermat(Xi,Zmesh,type,par,scale);
    %cond(A_xz)
    KXX=kermat(Xi,Xi,type,par,scale);
    

    condNum(nn) = cond(KXX);
%     L = chol(KXX);
%     l = 1./diag(L);
%     L = L*diag(l);
% 
%     condNew2(nn) = cond(L);

    %cond(A_xx)
    
    % Berechnung der u_r
    
    % ueber Inverse
%     u=kernel(distsq(X,X(end,:)));
%     Ainv=invUpdate(A,Ainv,u(1:end-1,:),u(end,:));
%     U1=Ainv*A_xz;
    
    % ueber LGS
    U=KXX\KXZ;
    ur(nn,:)=U(nn,:);
    

    index=[index,I];
    condInt(nn)=cond(ur(1:nn,index));
    
    % Ende Berechnung der u_r
    
    % bestimmen der Zeilensummennorm
    % der Matrix bei Newtonbasis
    newtonnorm(nn)=norm(inv(ur(1:nn,index)),inf);
    
    % bestimmen der Zeilensummennorm
    % der Inversen der Interpolationsmatrix
    intpolnorm(nn)=norm(inv(KXX),inf);
    
    % Berechnung der Powerfunktion P_r
    null=kermat(1,1,type,par,scale)*ones(len,1);
    P=null-reshape(sum(U.*KXZ,1),len,1);
    
    % Koeffizient lambda bestimmen
    lambda(nn)=yMesh(I)-lambda'*ur(:,I);
    
    % Interpolante bestimmen
    fn=lambda'*ur;
    err=fn'-yMesh;
    error(nn)=max(abs(err));
    norm_fs(nn) = sqrt(lambda'*lambda);
    
    % neuen Punkt hinzunehmen
    
    % ueber Powerfunktion
    [maxP(nn),I]=max(abs(P));
    if maxP(nn) < tol
        break;
    end
    % ueber Fehlerfunktion
    %[maxi,I]=max(abs(err));
    Xi=[Xi;Zmesh(I(1),:)];
end

if ~isempty(Xi(1:end-1,:))
    Xi=Xi(1:end-1,:);
end
norm_fs = max(norm_fs);

V = ur(1:nn,index)';
 % Punktmenge abspeichern mit Gridweite
 %save w2scale5peaks.mat X index h -ascii
 %save 'powerpoints1.mat' X -ascii
 %save 'powerpoints2.mat' index -ascii
 %save 'powerpoints3.mat' h -ascii
 %plot(X(1:end,1),X(1:end,2),'r+')
end