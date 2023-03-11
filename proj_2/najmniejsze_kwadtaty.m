clear endt start u y_p

dane_tr = load("dane.txt");
u_tr = dane_tr(:,1);
y_tr = dane_tr(:,2);

M = [];
for i=6:length(y_tr)
    row = [u_tr(i-4) u_tr(i-5) -y_tr(i-1) -y_tr(i-2)];
    M = [M;row];
end

y_tr = y_tr(6:end);  

b = M\y_tr;

rng(5)

endt = 3000;
start = 6;
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
y_p(1:endt) = 0;
Err = 0;

for k=start:endt
    g_1 = g1(u(k-4));
    x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
    x2(k) = -alpha2*x1(k-1)+betha2*g_1;
    y(k) = g2(x1(k));
    y_p(k) = [u(k-4),u(k-4-1),-y_p(k-1),-y_p(k-2)]*b;
    Err = Err + (y_p(k)-y(k))^2;
end

disp(Err)

figure(2)
hold on
stairs(y_p)
title("Wyjście weryfikujacy")
stairs(y)
legend("Przewidywane", "Rzeczywiste")
% print("wyjscie2.eps","-depsc","-r400")

figure(3)
stairs(u)
title("Sterowanie")
% print("sterowanie2.eps","-depsc","-r400")