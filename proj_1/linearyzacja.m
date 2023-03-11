clear all
close all

global Fd h1_LL h2_LL Fd_l F1_l C1 C2 alpha2 alpha1 T

C1=0.35;
C2 = 0.3;
alpha1 = 20;
alpha2 = 22;

T=5;
tau = round(160/T);



start = max(50, tau+2);
endt = 1000;

values = [42, 52, 62, 72, 82, 92];

for ff = values
F1(1:1:start+4) = 72;
F1(start+5:1:endt)=ff;
Fd = 15;


h1(1:1:start)=18.9225;
h1(start+1:1:endt)=0;
h2(1:1:start)=15.6402;
h2(start+1:1:endt)=0;

h1_l(1:1:start)=18.9225;
h1_l(start+1:1:endt)=0;
h2_l(1:1:start)=15.6402;
h2_l(start+1:1:endt)=0;

F1_l = 72;
Fd_l = 15;
h1_LL = 18.914;
h2_LL = 15.6248;


for k = start:1:endt

    [h1_prov, h2_prov] = func(h1(k-1), h2(k-1), F1(k-1-tau));
    [h1_mid, h2_mid] = func(h1(k-1)+h1_prov*T/2, h2(k-2)+h2_prov*T/2, F1(k-1-tau));
    h1(k) = h1(k-1) + T*h1_mid ;
    h2(k) = h2(k-1) + T*h2_mid;
    
    [h1l_prov, h2l_prov] = func_lin(h1_l(k-1), h2_l(k-1), F1(k-1-tau));
    [h1l_mid, h2l_mid] = func_lin(h1_l(k-1)+h1l_prov*T/2, h2_l(k-2)+h2l_prov*T/2, F1(k-1-tau));
    h1_l(k) = h1_l(k-1)+ h1l_mid*T;
    h2_l(k) = h2_l(k-1) + h2l_mid*T;
end
figure(1)
hold on
plot(h1)
   
figure(2)
hold on
plot(h2)
    
figure(3)
hold on
plot(h1_l)
    
figure(4)
hold on
plot(h2_l)
end

figure(1)
title("h1")
legend("F1 = "+flip(values))
print("h1.eps","-depsc","-r400")

figure(2)
title("h2")
legend("F1 = "+flip(values))
print("h2.eps","-depsc","-r400")

figure(3)
title("h1 zlinearyzowane")
legend("F1 = "+flip(values))
print("h1l.eps","-depsc","-r400")

figure(4)
title("h2 zlinearyzowane")
legend("F1 = "+flip(values))
print("h2l.eps","-depsc","-r400")

function [dh1, dh2] = func_lin(h1_l, h2_l, F1)
    global Fd h1_LL h2_LL Fd_l F1_l C1 C2 alpha2 alpha1
    el1 = (F1 + Fd - alpha1 * sqrt(h1_LL)) / (3*C1*h1_LL^2);
    el2 = ((h1_l - h1_LL) * (4*(F1_l + Fd_l) - 3*alpha1*sqrt(h1_LL))) / (6*C1*h1_LL^3);
    dh1 = (el1 - el2);
    el3 = (alpha1 * sqrt(h1_LL) - alpha2 * sqrt(h2_LL)) / (3*C2*h2_LL^2);
    el5 = (h1_l-h1_LL) * (alpha1/(6*sqrt(h1_LL)*h2_LL^2*C2));
    el4 = (h2_l-h2_LL) * (3*alpha2*sqrt(h2_LL) - 4*alpha1*sqrt(h1_LL)) / (6*C2*h2_LL^3);
    dh2 = (el3 +el5+ el4);
end

function [dh1, dh2] = func(h1, h2, F1)
    global C1 alpha1 C2 alpha2 Fd
    dh1 = (1/(3*C1*h1^2)) * (F1 + Fd - alpha1*sqrt(h1));
    dh2 =  (1/(3*C2*h2^2)) * (alpha1*sqrt(h1) - alpha2*sqrt(h2));
end

