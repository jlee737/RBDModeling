x1 = linspace(0, 0.01, 10000);
x2 = linspace(0,1,10000);
y1 = betapdf (x, 10, 1.3);
y2 = normpdf(x, 10^-5, 10^-5);


plot(x,y1)
hold on
plot(x,y2)
hold off
legend()

%betapdf(Kd, a, b)