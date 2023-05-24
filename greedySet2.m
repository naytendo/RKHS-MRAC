function [Xi,error,maxP,condNum,condInt,V] = greedySet2(tol,Omega,maxN,type,par,scale,fBound)

Z = Omega;
mu = 3;
yMesh = (1+mu)*Z(:,2)-mu*Z(:,1).^2.*Z(:,2);
len=size(Z,1);

% Startmenge [X,Y]
I=floor(len/2);
Xi=Z(I,:);
lambda=zeros(maxN,1);
ur=zeros(maxN,size(Z,1));
A=[];
Ainv=[];
error=zeros(maxN,1);
maxP=zeros(maxN,1);
condNum=zeros(maxN,1);
condInt=zeros(maxN,1);
intpolnorm=zeros(maxN,1);
norm_fs=zeros(maxN,1);
newtonnorm=zeros(maxN,1);
index=[];
beta=zeros(maxN); % Koeffizienten zur rekursiven Berechnung des OGS
%gamma=zeros(maxN,maxN);
ur2=zeros(maxN,size(Z,1));

for nn=1:maxN
    % Berechnung von A_xz und A_xx
    A_xz=kermat(Xi,Z,type,par,scale);
    %cond(A_xz)
    A_xx=kermat(Xi,Xi,type,par,scale);
    %cond(A_xx)
    condNum(nn) = cond(A_xx);
    % Berechnung der u_r
    
    % ueber Inverse
%     u=kernel(distsq(X,X(end,:)));
%     Ainv=invUpdate(A,Ainv,u(1:end-1,:),u(end,:));
%     U1=Ainv*A_xz;
    
    % ueber LGS
    U=A_xx\A_xz;
    ur(nn,:)=U(nn,:);
    
    % Ende Berechnung der u_r

    condInt(nn)=cond(ur(1:nn,index));
    % bestimmen der Zeilensummennorm
    % der Matrix bei Newtonbasis
    index=[index,I];
    newtonnorm(nn)=norm(inv(ur(1:nn,index)),inf);
    
    % bestimmen der Zeilensummennorm
    % der Inversen der Interpolationsmatrix
    intpolnorm(nn)=norm(inv(A_xx),inf);
    
    % Berechnung der Powerfunktion P_r
    null=kermat(1,1,type,par,scale)*ones(len,1);
    P=null-reshape(sum(U.*A_xz,1),len,1);
    
    % Koeffizient lambda bestimmen
    lambda(nn)=yMesh(I)-lambda'*ur(:,I);

    norm_fs(nn) = sqrt(lambda'*lambda);
    
    % Interpolante bestimmen
    fn=lambda'*ur;
    err=fn'-yMesh;
    error(nn)=max(abs(err));
    % neuen Punkt hinzunehmen
    
    % ueber Powerfunktion
    [maxPUn,I]=max(abs(P));
    maxP(nn) = maxPUn; % need to normalize
    if maxP(nn)*fBound < tol
        break;
    end
    % ueber Fehlerfunktion
    %[maxi,I]=max(abs(err));
    Xi=[Xi;Z(I(1),:)];
end
if size(Xi,1) > 1
    Xi=Xi(1:end-1,:);
end
V = real(ur(1:nn,index)');
