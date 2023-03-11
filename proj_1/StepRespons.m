function respons = StepRespons(skok_sterowania, h1l, h2l, rysowanie)
%     skok_sterowania = 73;
%     h1l = 18.914;
%     h2l = 15.6248;
%     rysowanie = 1;
    
    start = 200;
    endt = 2000;
    
    global Fd h1_LL h2_LL Fd_l F1_l C1 C2 alpha2 alpha1 T

    T=5;
    tau = round(160/T);
    start = tau+2;

    C1 = 0.35;
    C2 = 0.3;
    alpha1 = 20;
    alpha2 = 22;
    F1(1:1:start) = 72;
    F1(start:1:endt) = skok_sterowania;
    Fd = 15;

    h1(1:1:start)=18.9225;
    h1(start+1:1:endt)=0;
    h2(1:1:start)=15.6402;
    h2(start+1:1:endt)=0;

    h1_l(1:1:start)=18.9225;
    h1_l(start+1:1:endt)=0;
    h2_l(1:1:start)=15.6402;
    h2_l(start+1:1:endt)=0;
    


    h1_LL = h1l;
    h2_LL = h2l;
    F1_l = skok_sterowania;
    Fd_l = Fd;

    for k = start:1:endt
        
        [h1l_prov, h2l_prov] = func_lin(h1_l(k-1), h2_l(k-1), F1(k-1-tau));
        [h1l_mid, h2l_mid] = func_lin(h1_l(k-1)+h1l_prov*T/2, h2_l(k-2)+h2l_prov*T/2, F1(k-1-tau));
        h1_l(k) = h1_l(k-1)+ h1l_mid*T;
        h2_l(k) = h2_l(k-1) + h2l_mid*T;
        h1_l(k) = max(0, h1_l(k));
        h2_l(k) = max(0, h2_l(k));

    end
    
    skok = [];
    i = 1;
    for k=1:1:length(h2_l)
        skok(i)= (h2_l(k)-15.6402)/abs((F1(start)-F1(start-1))); % bo sterowanie w procentach
        i = i + 1;
    end
    respons = skok(start:1:endt);
    if rysowanie == 1
        stairs(respons);
        title("Przeskalowana odpowiedź zlinearyzowanego obiektu")
        print("odp_skok_przeskalowana.eps","-depsc","-r400")
    end
 end


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