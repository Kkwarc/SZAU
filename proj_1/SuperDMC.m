function [E, h1, h2, h2_zad, F1, Fd]=SuperDMC(wektor, zaklocenia, rysowanie)

    T=5;
    tau = round(160/T);
    E=0;
    
    s_load = StepRespons(73, 18.9225, 15.6402, 0); % poj. odp. skokowa
    
    
    D = 500;
    s = s_load(1:1:D);

    
    start= D+1;
    endt=5000; %koniec symulacji
    poczatek = start; %chwila k w której zmienia sie wartość zadana

    
    %Definicja horyzontów i parametrów
    N = round(wektor(1));
    N_u = round(wektor(2));
    lambda = round(wektor(3));

    
    u_min = 0.5*72;
    u_max = 1.5*72;

    
    %deklaracja potrzebnych wektorów
    h1(1:1:start)=18.9225;
    h1(start+1:1:endt)=0;
    h2(1:1:start)=15.6402;
    h2(start+1:1:endt)=0;

    C1=0.35;
    C2 = 0.3;
    alpha1 = 20;
    alpha2 = 22;

    

    F1(1:1:start+1) = 72;
    F1(start+2:1:endt) = 0;

    %Obliczenie części macierzy DMC
    M = zeros(N, N_u);
    for column=1:N_u
        for row=1:N
           if (row>=column)
             if(row-column+1<=D)
                    M(row,column)=s(row-column+1);
                  else
                   M(row,column)=s(D);
             end
          end
        end
    end

    K = (M'*M+lambda*eye(N_u, N_u))^(-1)*M';

    M_p = zeros(N, D-1);
    for column=1:(D-1)
        for row=1:N
            if row + column > D
                if column>D
                    M_p(row, column) = 0;
                else
                    M_p(row, column) = s(D) - s(column);
                end
            else
                M_p(row, column) = s(row + column) - s(column);
            end
        end
    end

    
    %skok wartości zadanej
    h2_zad(1:poczatek)=15.6402; 
    h2_zad(poczatek+1:1500)=31;
    h2_zad(1501:2500)=11;
    h2_zad(2501:3500)=25;
    h2_zad(3501:endt)= 8;
    
    if zaklocenia == 0
        % bez zakłóceń
        Fd(1:1:endt) = 15;
    elseif zaklocenia == 1
        % z zakłoceniami
        Fd(1:1000)=15;
        Fd(1001:2000)=24;
        Fd(2001:4000)=8;
        Fd(4001:endt)=0;
    end

    
    %inicjalizacja pozostałych potrzebnych macierzy
    DU_p = zeros(D-1, 1);

    for k=start:endt
        %symulacja obiektu
        h1(k) = h1(k-1) + T * (1/(3*C1*h1(k-1)^2)) * (F1(k-1-tau) + Fd(k) - alpha1*sqrt(h1(k-1)));
        if h1(k)<0
            h1(k) = 0;
        end
        h2(k) = h2(k-1) + T * (1/(3*C2*h2(k-1)^2)) * (alpha1*sqrt(h1(k-1)) - alpha2*sqrt(h2(k-1)));
        if h2(k)<0
            h2(k)=0;
        end

        %Obliczenie DU_p
        for d=1:(D-1)
            DU_p(d) = F1(k-d) - F1(k-d-1);
        end

        %Pomiar wyjścia
        Y = ones(N, 1) * h2(k);

        %Obliczenie Y_0
        yo = M_p * DU_p + Y;

        Y_zad = ones(N, 1) * h2_zad(k);

        %Obliczenie sterowania
        DU = K * (Y_zad - yo);
        F1(k) = F1(k-1) + DU(1);
        F1(k) = min(F1(k), u_max);
        F1(k) = max(F1(k), u_min);
        E = E + (h2_zad(k)-h2(k))^2; 
    end
    if rysowanie == 1
        subplot(2,1,1)
        stairs(h2)
        hold on
        stairs(h2_zad)
        title("Wyjście")
        legend("Wyjście", "Wartość zadana")
        xlabel("chwila k")
        ylabel("wartość wyjścia")
        hold off
        subplot(2,1,2)
        hold on
        stairs(F1)
        stairs(Fd)
        legend("sterowanie", "zakłócenia")
        xlabel("chwila k")
        ylabel("wartość sterowania")
        title("Sterowanie")
        print("DMC_single.eps","-depsc","-r400")
    end
end