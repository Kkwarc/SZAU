clear all
close all

trainf = 'trainlm';  % Levenberg-Marquardt
% trainf = 'trainscg';  % Scaled Conjugate Gradient

net = feedforwardnet(9, trainf);

endt = 3000;
start = 6;
alpha1 = -1.599028;
alpha2 = 0.632337;
betha1 = 0.010754;
betha2 = 0.009231;

u_min = -1;
u_max = 1;

rng(1)

ut(1:endt) = losowe_sterowanie(50, endt, -1, 1);
x1t(1:endt) = 0;
x2t(1:endt) = 0;
yt(1:endt) = 0;
for k=start:endt
    g_1t = g1(ut(k-4));
    x1t(k) = -alpha1*x1t(k-1)+x2t(k-1)+betha1*g_1t;
    x2t(k) = -alpha2*x1t(k-1)+betha2*g_1t;
    yt(k) = g2(x1t(k));
end

u4 = ut(2:end-4);
u5 = ut(1:end-5);
y1 = yt(5:end-1);
y2 = yt(4:end-2);
yl = yt(6:end);
dl = [u4; u5; y1; y2];

rng(5) % do zamiany przy probach uczenia
net = train(net,dl,yl);
disp("trained")

rng(5)

u(1:endt) = losowe_sterowanie(50, endt, -1, 1);
x1(1:endt) = 0;
x2(1:endt) = 0;
y(1:endt) = 0;
y_p_arx(1:endt) = 0;
y_p_oe(1:endt) = 0;
y_u_arx(1:endt) = 0;
y_u_oe(1:endt) = 0;
ErrTArx = 0;
ErrTOe = 0;
ErrUArx = 0;
ErrUOe = 0;

for k=start:endt
    g_1 = g1(u(k-4));
    x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
    x2(k) = -alpha2*x1(k-1)+betha2*g_1;
    y(k) = g2(x1(k));
    
    % ARX
    y_p_arx(k) = sim(net, [u(k-4); u(k-5); y(k-1); y(k-2)]);
    ErrTArx = ErrTArx + (y_p_arx(k)-y(k))^2;
    
    % OE
    y_p_oe(k) = sim(net, [u(k-4); u(k-5); y_p_oe(k-1); y_p_oe(k-2)]);
    ErrTOe = ErrTOe + (y_p_oe(k)-y(k))^2;
    
    % test dla danych uczących ARX
    y_u_arx(k) = sim(net, [ut(k-4); ut(k-5); yt(k-1); yt(k-2)]);
    ErrUArx = ErrUArx + (y_u_arx(k)-yt(k))^2;
    
    % test dla danych uczących OE
    y_u_oe(k) = sim(net, [ut(k-4); ut(k-5); y_u_oe(k-1); y_u_oe(k-2)]);
    ErrUOe = ErrUOe + (y_u_oe(k)-yt(k))^2;
end

figure(1)
subplot(2,1,1)
hold on
stairs(y_p_arx)
title(['Wyjście ARX - dane testujące, błąd: ', num2str(ErrTArx)])
stairs(y)
legend("Przewidywane", "Rzeczywiste")
subplot(2,1,2)
stairs(u)
title("Sterowanie")
print("Zad3ARX5.eps","-depsc","-r400")

figure(2)
subplot(2,1,1)
hold on
stairs(y_p_oe)
title(['Wyjście OE - dane testujące, błąd: ', num2str(ErrTOe)])
stairs(y)
legend("Przewidywane", "Rzeczywiste")
subplot(2,1,2)
stairs(u)
title("Sterowanie")
print("Zad3Oe5.eps","-depsc","-r400")

figure(3)
subplot(2,1,1)
hold on
stairs(y_u_arx)
title(['Wyjście ARX- dane uczące, błąd: ', num2str(ErrUArx)])
stairs(yt)
legend("Przewidywane", "Rzeczywiste")
subplot(2,1,2)
stairs(ut)
title("Sterowanie")
print("Zad3Ucz5ARX.eps","-depsc","-r400")

figure(4)
subplot(2,1,1)
hold on
stairs(y_u_oe)
title(['Wyjście OE - dane uczące, błąd: ', num2str(ErrUOe)])
stairs(yt)
legend("Przewidywane", "Rzeczywiste")
subplot(2,1,2)
stairs(ut)
title("Sterowanie")
print("Zad3Ucz5OE.eps","-depsc","-r400")