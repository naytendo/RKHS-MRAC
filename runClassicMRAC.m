close all; clear; clc;

a_plant=-10; b_plant=+10;
ar=-20; br=20;
ga=10; gb=10; gf=10;

x0=1; xr0=0;
alpha0=0; beta0=0;sigma=0.1; N=3;

basiscenter=linspace(-10,10,N);

knlfun=@(x,bsc) exp(-(x-bsc)^2/(2*(1/N)^2));
k_xi_xi=zeros(length(basiscenter));
for i=1:length(basiscenter)
    for j=1:length(basiscenter)
        k_xi_xi(i,j)=knlfun(basiscenter(i),basiscenter(j));
    end
end

open classicMRAC.slx
sim('classicMRAC.slx')