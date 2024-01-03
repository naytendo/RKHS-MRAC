
x1 = linspace(-2,-1,100);
x2 = linspace(-1,1,100);
x3 = linspace(1,2,100);
y1 = -1*ones(length(x1));
y3 = 1*ones(length(x1));
y2 = x2;

figure()
plot(x1,y1,'k','linewidth',2)
hold on
plot(x2,y2,'k','linewidth',2)
plot(x3,y3,'k','linewidth',2)
grid on
set(gca,'fontSize',20)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
ylim([-1.5 1.5])
line([-1 -1],[-1.5 1.5])
line([1 1],[-1.5 1.5])