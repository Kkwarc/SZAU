clear all

C1=0.35;
C2 = 0.3;
alpha1 = 20;
alpha2 = 22;
F1 = 72;
Fd = 15;

dh = @(t,h) [(1/(3*C1*h(1)^2))*(F1 + Fd - alpha1*sqrt(h(1))); (1/(3*C2*h(2)^2))*(alpha1*sqrt(h(1))-alpha2*sqrt(h(2)))];

Y = ode45(dh, [0, 1000], [2, 3]);

plot(Y.x, Y.y)
title("Punkt pracy")
% legend("h1", "h2", "Location", "SE")
print("punkt_pracy.eps","-depsc","-r400")