%PID------------------------------------------------------------------
function [E, Y, y_zad, U]=Pid(wektor)
    T = 1; % okres próbkowania 
    K=wektor(1); % człon proporcjonalny
    Ti=wektor(2); % człon całkujący
    Td=wektor(3); % człon różniczkujący 
%     disp(wektor)

    % parametry dyskretnego regulatora PID
    r2 = (K*Td)/T;
    r1 = K*(T/(2*Ti) - 2*(Td/T) - 1);
    r0 = K*(1+ T/(2*Ti) + Td/T);
    Ts=5000+12;
    
    %warunki początkowe
    alpha1 = -1.599028;
    alpha2 = 0.632337;
    betha1 = 0.010754;
    betha2 = 0.009231;

    u_min = -1;
    u_max = 1;

    U(1:Ts)=0; Y(1:Ts)=0; e(1:Ts)=0; x1(1:Ts)=0; x2(1:Ts)=0;
    y_zad(1:12)=0; 
    y_zad(12+1:1500+12)=2;
    y_zad(1501+12:12+2500)=-1;
    y_zad(2501+12:12+3500)=3;
    y_zad(3501+12:Ts)= 0;
    E=0;
    
    for k=12:Ts %główna pętla symulacyjna
        %symulacja obiektu
        g_1 = g1(U(k-4));
        x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
        x2(k) = -alpha2*x1(k-1)+betha2*g_1;
        Y(k) = g2(x1(k));
        %uchyb regulacji
        e(k)=y_zad(k)-Y(k);
        %sygnał sterujący regulatora PID
        U(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+U(k-1);

        if U(k)>u_max
            U(k)=u_max;
        end
        if U(k)<u_min
            U(k)=u_min;
        end   
        E = E + (y_zad(k)-Y(k))^2;
    end
end