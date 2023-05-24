x0 = [2;3];
tspan = [0 100];

options = odeset('RelTol',1e-7,'AbsTol',1e-8);
x = ode45(@vanDerPol,tspan,x0,options);
x = x.y;
plot(x(1,1),x(2,1),'ko','linewidth',1.5);
hold on
plot(x(1,:),x(2,:),'k','linewidth',1.5);

xlabel('$x_1$','interpreter','latex','fontsize',20)
ylabel('$x_2$','interpreter','latex','fontsize',20)
grid on
legend('$x_0$','$x(t)$','interpreter','latex','fontsize',20)