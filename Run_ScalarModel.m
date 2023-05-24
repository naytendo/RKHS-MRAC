close all; clear; clc;

astar=-10; bstar=+10;
ar=-20; br=20;
ga=10; gb=10; gf=10;

x0=0; xr0=0;
alpha0=0; beta0=0;sigma=0.1; N=12;
kx0=0; kr0=0;


%% RKHS Preparation
basiscenter=linspace(-10,10,N)';
hype=1/(2*2^2); l=2; sig=2;

% Exponential Kernel
knlfune=@(x,bsc) exp(-hype*(x-bsc).^2);
% Matern 3/2 Kernel
knlfun3m=@(x,bsc) sig^2.*(1+sqrt(3)./l.*abs(x-bsc)).*exp(-sqrt(3).*abs(x-bsc)./l);
% Matern 5/2 Kernel
knlfun5m=@(x,bsc) sig^2.*(1+sqrt(5)/l.*abs(x-bsc)+5/(3*l^2).*abs(x-bsc).^2).*exp(-sqrt(5).*abs(x-bsc)./l);
k_xi_xie=zeros(N);
k_xi_xi3m=zeros(N);
k_xi_xi5m=zeros(N);
for i=1:N
    for j=1:N
        k_xi_xie(i,j)=knlfune(basiscenter(i),basiscenter(j));
    end
end
condnume=cond(k_xi_xie);

inv_k_xi_xie=inv(k_xi_xie);

Thetastare=inv_k_xi_xie*(10*tanh(basiscenter));

for i=1:N
    for j=1:N
        k_xi_xi3m(i,j)=knlfun3m(basiscenter(i),basiscenter(j));
    end
end
condnum3m=cond(k_xi_xi3m);

inv_k_xi_xi3m=inv(k_xi_xi3m);

Thetastar3m=inv_k_xi_xi3m*(10*tanh(basiscenter));

for i=1:N
    for j=1:N
        k_xi_xi5m(i,j)=knlfun5m(basiscenter(i),basiscenter(j));
    end
end
condnum5m=cond(k_xi_xi5m);

inv_k_xi_xi5m=inv(k_xi_xi5m);

Thetastar5m=inv_k_xi_xi5m*(10*tanh(basiscenter));


open ScalarModel.slx
sim('ScalarModel.slx')