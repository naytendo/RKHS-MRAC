clear;
close all;

global cscale;
global kertype;

cscale=5;
kertype='w2'; % Basisfunktion
%kertype='g'; % Basisfunktion
fkt='funcR2'; % Datenfunktion

maxN=80; % maximale Anzahl Iterationen
h=0.05; % Gitterweite, feines Gitter
[Z1,Z2]=meshgrid(-3:h:3,-3:h:3);
Z=[Z1(:),Z2(:)];
f=evalFuncRd(fkt,Z');
len=size(Z,1);

% Startmenge [X,Y]
I=floor(len/2);
X=Z(I,:);
lambda=zeros(maxN,1);
ur=zeros(maxN,size(Z,1));
A=[];
Ainv=[];
error=zeros(maxN,1);\
norm_fs=zeros(maxN,1);
intpolnorm=zeros(maxN,1);
newtonnorm=zeros(maxN,1);
index=[];
beta=zeros(maxN); % Koeffizienten zur rekursiven Berechnung des OGS
%gamma=zeros(maxN,maxN);
ur2=zeros(maxN,size(Z,1));

for r=1:maxN
    r
    % Berechnung von A_xz und A_xx
    A_xz=kernel(distsq(X,Z));
    %cond(A_xz)
    A_xx=kernel(distsq(X,X));
    %cond(A_xx)
    
    % Berechnung der u_r
    
    % ueber Inverse
%     u=kernel(distsq(X,X(end,:)));
%     Ainv=invUpdate(A,Ainv,u(1:end-1,:),u(end,:));
%     U1=Ainv*A_xz;
    
    % ueber LGS
    U=A_xx\A_xz;
    ur(r,:)=U(r,:);
    
    % Ende Berechnung der u_r
    
    % bestimmen der Zeilensummennorm
    % der Matrix bei Newtonbasis
    index=[index,I];
    newtonnorm(r)=norm(inv(ur(1:r,index)),inf);
    
    % bestimmen der Zeilensummennorm
    % der Inversen der Interpolationsmatrix
    intpolnorm(r)=norm(inv(A_xx),inf);
    
    % Berechnung der Powerfunktion P_r
    null=kernel(0)*ones(len,1);
    P=null-reshape(sum(U.*A_xz,1),len,1);
    
    % Koeffizient lambda bestimmen
    lambda(r)=evalFuncRd(fkt,Z(I,:)')-lambda'*ur(:,I);
    
    % Interpolante bestimmen
    fn=lambda'*ur;
    err=fn'-f;
    error(r)=max(abs(err));
    norm_fs(r) = sqrt(lambda'*lambda);
    
    % neuen Punkt hinzunehmen
    
    % ueber Powerfunktion
    [maxi,I]=max(abs(P));
    % ueber Fehlerfunktion
    %[maxi,I]=max(abs(err));
    X=[X;Z(I(1),:)];
end
X=X(1:end-1,:);
norm_fs = max(norm_fs);
 % Punktmenge abspeichern mit Gridweite
 %save w2scale5peaks.mat X index h -ascii
 %save 'powerpoints1.mat' X -ascii
 %save 'powerpoints2.mat' index -ascii
 %save 'powerpoints3.mat' h -ascii
 save powerpoints.mat X index h
 %plot(X(1:end,1),X(1:end,2),'r+')