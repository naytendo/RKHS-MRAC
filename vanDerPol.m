function dxdt = vanDerPol(t,x)
mu = 3;

dxdt = [x(2); mu*(1-x(1)^2)*x(2)-x(1)];

end