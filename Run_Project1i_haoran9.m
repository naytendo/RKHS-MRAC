close all; clear; clc;

astar=-10; bstar=+10;
ar=-20; br=20;
gx=10; gr=10;sigma=0.1;

x0=1; xr0=0;
kx0=0; kr0=0;

kxstar=(ar-astar)/bstar; krstar=br/bstar;

open Project1i_haoran9.slx
sim('Project1i_haoran9.slx')