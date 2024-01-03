x = linspace(0.1,100,1000);

y = exp(-0.04*x).*cos(1/20*x)+exp(-0.0005*x)+0.2*cos(.8*x).*exp(-0.1*x);

figure()
plot(x,y,'k','linewidth',2)
grid on
set(gca,'fontSize',20)
set(gca,'YTickLabel',[])


line([11.4 17],[1.34 1.34])
line([17 20],[1.34 1.34])
line([13.4 13.4],[1.34 0])
line([17 17],[1.34 0])
line([22 22],[1.34 0])
line([22 27],[1.048 1.048])
line([27 30],[1.048 1.048])
line([27 27],[1.048 0])
line([32 32],[1.048 0])
line([32 90],[0.86 0.86])
text(8,1,'$t_1$','interpreter','Latex','fontsize',25)
text(8,1,'$t_2$','interpreter','Latex','fontsize',25)
set(gca,'XTickLabel',[])

text(8,1,'$t_1$','interpreter','Latex','fontsize',25)
text(8,1.2,'$t_2$','interpreter','Latex','fontsize',25)

text(8,1.4,'$t_3$','interpreter','Latex','fontsize',25)