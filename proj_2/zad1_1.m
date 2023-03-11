clear all
close all

endt = 200;
start = 5;
alpha1 = -1.599028;
alpha2 = 0.632337;
betha1 = 0.010754;
betha2 = 0.009231;

u_min = -1;
u_max = 1;

U(1:(u_max-u_min)/0.1+1) = [u_min:0.1:u_max];
Y(1:(u_max-u_min)/0.1+1) = 0;

for i=1:length(U)
    u(1:start) = 0;
    u(start:endt) = U(i);
    x1(1:endt) = 0;
    x2(1:endt) = 0;
    y(1:endt) = 0;
    for k=start:endt
        g_1 = g1(u(k-4));
        x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
        x2(k) = -alpha2*x1(k-1)+betha2*g_1;
        y(k) = g2(x1(k));
    end
    Y(i) = y(endt);
end

figure(1)
plot([u_min:0.1:u_max], Y)
title("Charakterystyka statyczna")
print("CharakterystykaStatyczna.eps","-depsc","-r400")