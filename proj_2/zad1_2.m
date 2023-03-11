clear all
close all

rng(1)

endt = 3000;
start = 5;
alpha1 = -1.599028;
alpha2 = 0.632337;
betha1 = 0.010754;
betha2 = 0.009231;

u_min = -1;
u_max = 1;

u(1:endt) = losowe_sterowanie(100, endt, -1, 1);
x1(1:endt) = 0;
x2(1:endt) = 0;
y(1:endt) = 0;
for k=start:endt
    g_1 = g1(u(k-4));
    x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
    x2(k) = -alpha2*x1(k-1)+betha2*g_1;
    y(k) = g2(x1(k));
end

figure(2)
stairs(y)
title("Wyjście")
figure(3)
stairs(u)
title("Sterowanie")
print("LosowePobudzenie.eps","-depsc","-r400")