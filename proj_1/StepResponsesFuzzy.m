function responses = StepResponsesFuzzy(liczba_regulatorow, bounds, rysowanie)
    responses = cell(1, liczba_regulatorow);
    bound_min = bounds(1);
    bound_max = bounds(2);

    Y = (bound_min+(bound_max-bound_min)/(liczba_regulatorow+1)):((bound_max-bound_min)/(liczba_regulatorow+1)):(bound_max-(bound_max-bound_min)/(liczba_regulatorow+1));
    i = 1;
    global Fd h1_LL h2_LL Fd_l F1_l C1 C2 alpha2 alpha1 T
    for y=Y

        C1=0.35;
        C2 = 0.3;
        alpha1 = 20;
        alpha2 = 22;

        T=5;
        tau = round(160/T);
        h2_LL = y;
        Fd_l = 15;
        F1_l = sqrt(h2_LL)*alpha2-Fd_l;
        h1_LL = ((F1_l+Fd_l)/alpha1)^2;
  
        start = max(2, tau+3);
        endt = start+100*50/T;
        
        F1(1:1:start) = F1_l;
        F1(start:1:endt)=F1_l+1;
        Fd = 15;
         
        h1(1:1:start)=h1_LL;
        h1(start+1:1:endt)=0;
        h2(1:1:start)=h2_LL;
        h2(start+1:1:endt)=0;
        
        h1_l(1:1:start)=h1_LL;
        h1_l(start+1:1:endt)=0;
        h2_l(1:1:start)=h2_LL;
        h2_l(start+1:1:endt)=0;

        for k = start:1:endt
            [h1l_prov, h2l_prov] = func_lin(h1_l(k-1), h2_l(k-1), F1(k-1-tau));
            [h1l_mid, h2l_mid] = func_lin(h1_l(k-1)+h1l_prov*T/2, h2_l(k-2)+h2l_prov*T/2, F1(k-1-tau));
            h1_l(k) = h1_l(k-1)+ h1l_mid*T;
            h2_l(k) = h2_l(k-1) + h2l_mid*T;
            h1_l(k) = max(0, h1_l(k));
            h2_l(k) = max(0,h2_l(k));
        end
        
        stat = (h2_l(start+1:1:endt)-h2_LL)/abs(F1(1)-F1(start));
        
        responses{i} = stat;
        i=i+1;
        if rysowanie == 1
            figure(8)
            hold on
            stairs(stat)
        end
    end
    if rysowanie == 1
        figure(8)
        legend("Odpowiedź skokowa dla modelu "+(1:1:i), "Location", "SouthEast")
        title("Opowiedzi skokowe dla modeli lokalnych w ilości: "+liczba_regulatorow)
        print("model_rozmyty_5.eps","-depsc","-r400")
    end
end

function [h1, h2] = func_lin(h1_l, h2_l, F1)
    global Fd h1_LL h2_LL Fd_l F1_l C1 C2 alpha2 alpha1
    el1 = (F1 + Fd - alpha1 * sqrt(h1_LL)) / (3*C1*h1_LL^2);
    el2 = ((h1_l - h1_LL) * (4*(F1_l + Fd_l) - 3*alpha1*sqrt(h1_LL))) / (6*C1*h1_LL^3);
    h1 = (el1 - el2);
    el3 = (alpha1 * sqrt(h1_LL) - alpha2 * sqrt(h2_LL)) / (3*C2*h2_LL^2);
    el5 = (h1_l-h1_LL) * (alpha1/(6*sqrt(h1_LL)*h2_LL^2*C2));
    el4 = (h2_l-h2_LL) * (3*alpha2*sqrt(h2_LL) - 4*alpha1*sqrt(h1_LL)) / (6*C2*h2_LL^3);
    h2 = (el3 +el5+ el4);
end