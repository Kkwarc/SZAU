clear all
global Fd h1_LL h2_LL Fd_l F1_l C1 C2 alpha2 alpha1 T
liczba_regulatorow = 2;


Fd_l = 15;  
Fd = 15;

C1=0.35;
C2 = 0.3;
alpha1 = 20;
alpha2 = 22;

T=5;
tau = round(160/T);

bound_min = 5;
bound_max = 35;

start = max(4, tau+4);
endt = 500;

ef1 = zeros(1, liczba_regulatorow);
ef2 = zeros(1, liczba_regulatorow);
ef3 = zeros(1, liczba_regulatorow);
ef4 = zeros(1, liczba_regulatorow);
ef5 = zeros(1, liczba_regulatorow);

h2_FF = (bound_min+(bound_max-bound_min)/(liczba_regulatorow+1)):((bound_max-bound_min)/(liczba_regulatorow+1)):(bound_max-(bound_max-bound_min)/(liczba_regulatorow+1));

F1_f = alpha2*sqrt(h2_FF)-Fd_l;

h1_FF = ((F1_f+Fd_l)./alpha1).^2;

width = bound_max-bound_min;
functions = MembershipFunction(liczba_regulatorow, 'tr', bound_min, bound_max, 0);

ff=90;
F1(1:1:start+4) = 50;
F1(start+5:1:endt)=ff;

h1_f(1:1:start)=10.56;
h1_f(start+1:1:endt)=0;
h2_f(1:1:start)=8.73;
h2_f(start+1:1:endt)=0;

h1(1:1:start)=10.56;
h1(start+1:1:endt)=0;
h2(1:1:start)=8.73;
h2(start+1:1:endt)=0;

h1_l(1:1:start)=10.56;
h1_l(start+1:1:endt)=0;
h2_l(1:1:start)=8.73;
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


    for i=1:1:liczba_regulatorow
        ef1(i) = (F1(k-1-tau) + Fd - alpha1 * sqrt(h1_FF(i))) / (3*C1*h1_FF(i)^2);
        ef2(i) = ((h1_f(k-1) - h1_FF(i)) * (4*(F1_f(i) + Fd_l) - 3*alpha1*sqrt(h1_FF(i)))) / (6*C1*h1_FF(i)^3);
        ef3(i) = (alpha1 * sqrt(h1_FF(i)) - alpha2 * sqrt(h2_FF(i))) / (3*C2*h2_FF(i)^2);
        ef5(i) = (h1_f(k-1)-h1_FF(i)) * (alpha1/(6*sqrt(h1_FF(i))*h2_FF(i)^2*C2));
        ef4(i) = (h2_f(k-1)-h2_FF(i)) * (3*alpha2*sqrt(h2_FF(i)) - 4*alpha1*sqrt(h1_FF(i))) / (6*C2*h2_FF(i)^3);
    end
    
    sum_mi = 0;
    for i=1:liczba_regulatorow
        sum_mi = sum_mi + functions{i}(h2_f(k-1));
    end
    
    EF1=0;
    EF2=0;
    EF3=0;
    EF4=0;
    EF5=0;
    
    for i=1:1:liczba_regulatorow
        EF1 = EF1 +functions{i}(h2_f(k-1))/sum_mi*ef1(i);
        EF2 = EF2 +functions{i}(h2_f(k-1))/sum_mi*ef2(i);
        EF3 = EF3 +functions{i}(h2_f(k-1))/sum_mi*ef3(i);
        EF4 = EF4 +functions{i}(h2_f(k-1))/sum_mi*ef4(i);
        EF5 = EF5 +functions{i}(h2_f(k-1))/sum_mi*ef5(i);
    end
    
    h1_f(k) = h1_f(k-1) + T * (EF1 - EF2);
    h2_f(k) = h2_f(k-1) + T * (EF3 +EF5+ EF4);
end

set(0, 'DefaultLineLineWidth', 1);
figure(1)
hold on
stairs(h2_f, 'r-');
stairs(h2, 'b-');
stairs(h2_l,'k-.');
title('2 modele lokalne - funkcje trójkątne')
legend('model rozmyty', 'model nieliniowy', 'model zlinearyzowany','Location','SE')
% print("model_rozmyty_2lok_gaus.eps","-depsc","-r400")

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
