function [E, h1, h2, h2_zad, F1, Fd]=FuzzyDmc(wektor, liczba_regulatorow, typ_funkcji, zaklocenia, rysowanie)
%DMC------------------------------------------------------------------

E=0;
T=5;
tau = round(160/T);


odp_skok = StepResponsesFuzzy(liczba_regulatorow, [5, 35], 0);
M = cell(1,liczba_regulatorow);


y_min = 5;
y_max = 35;


functions = MembershipFunction(liczba_regulatorow, typ_funkcji, y_min, y_max, 0);


D = 500;
start= D+1;
endt=5000+D; %koniec symulacji
poczatek = start; %chwila k w której zmienia sie wartość zadana zalecana > 80


%Definicja horyzontów i parametrów
N = wektor(1);
N_u = wektor(2);
lambda = wektor(3);


%Obliczenie części macierzy DMC
for i=1:liczba_regulatorow
    M{i} = zeros(N, N_u);
    for column=1:N_u
        for row=1:N
           if (row>=column)
             if(row-column+1<=D)
                    M{i}(row,column)=odp_skok{i}(row-column+1);
                  else
                   M{i}(row,column)=odp_skok{i}(D);
             end
          end
        end
    end
end

K = cell(1,liczba_regulatorow);
for i=1:liczba_regulatorow
    K{i} = (M{i}'*M{i}+lambda*eye(N_u, N_u))^(-1)*M{i}';
end


%deklaracja potrzebnych wektorów
h1(1:1:start)=18.914;
h1(start+1:1:endt)=0;
h2(1:1:start)=15.6548;
h2(start+1:1:endt)=0;

C1=0.35;
C2 = 0.3;
alpha1 = 20;
alpha2 = 22;

u_min = sqrt(y_min)*alpha2-15;
u_max = sqrt(y_max)*alpha2-15;

F1(1:1:start+1) = 72;
F1(start+2:1:endt) = 0;


M_p = cell(1, liczba_regulatorow);
% M{i} = zeros(N, N_u);
for i=1:liczba_regulatorow
    M_p{i} = zeros(N, D-1);
    for column=1:(D-1)
        for row=1:N
            if row + column > D
                if column>D
                    M_p{i}(row, column) = 0;
                else
                    M_p{i}(row, column) = odp_skok{i}(D) - odp_skok{i}(column);
                end
            else
                M_p{i}(row, column) = odp_skok{i}(row + column) - odp_skok{i}(column);
            end
        end
    end
end

    %skok wartości zadanej
    h2_zad(1:poczatek)=15.6548; 
    h2_zad(poczatek+1:1500+D)=31;
    h2_zad(1501+D:2500+D)=11;
    h2_zad(2501+D:3500+D)=25;
    h2_zad(3501+D:endt)= 8;

    if zaklocenia == 0
        % bez zakłóceń
        Fd(1:1:endt) = 15;
    elseif zaklocenia == 1
        % z zakłoceniami
        Fd(1:1000+D)=15;
        Fd(1001+D:2000+D)=24;
        Fd(2001+D:4000+D)=8;
        Fd(4001+D:endt)=0;
    end


%inicjalizacja pozostałych potrzebnych macierzy
DU_p = zeros(D-1, 1);


for k=start:endt
    %symulacja obiektu
    h1(k) = h1(k-1) + T * (1/(3*C1*h1(k-1)^2)) * (F1(k-1-tau) + Fd(k) - alpha1*sqrt(h1(k-1)));
    h1(k) = max(0,h1(k));
    
    h2(k) = h2(k-1) + T * (1/(3*C2*h2(k-1)^2)) * (alpha1*sqrt(h1(k-1)) - alpha2*sqrt(h2(k-1)));
    h2(k) = max(0,h2(k));

    %Obliczenie DU_p
    for d=1:(D-1)
        DU_p(d) = F1(k-d) - F1(k-d-1);
    end

    %Obliczenie sterowania
    sum_mi = 0;
    for i=1:liczba_regulatorow
        sum_mi = sum_mi + functions{i}(h2(k));
    end

    DU = cell(1,liczba_regulatorow);
    for i=1:liczba_regulatorow
         %Pomiar wyjścia
        Y = ones(N, 1) * h2(k);

        Y_zad = ones(N, 1) * h2_zad(k);
        
        DU{i} = (K{i} * (Y_zad -  M_p{i} * DU_p - Y));
    end
    
    u=0;
    for i =1:liczba_regulatorow
        u = u + (functions{i}(h2(k))/sum_mi)*DU{i}(1); 
    end

    F1(k) = F1(k-1) + u;
    F1(k) = max(F1(k), u_min+0.0001);
    F1(k) = min(F1(k), u_max-0.0001);
    
    E=E+(h2_zad(k)-h2(k))^2;
end
    if rysowanie == 1
        subplot(2,1,1)
        stairs(h2)
        hold on;
        stairs(h2_zad)
        title("Wyjście")
        legend("Wyjście", "Wartość zadana")
        xlabel("chwila k")
        ylabel("wartość wyjścia")
        hold off;

        subplot(2,1,2)
        hold on;
        stairs(F1)
        stairs(Fd)
        legend("sterowanie", "zakłócenia")
        xlabel("chwila k")
        ylabel("wartość sterowania")
        title("Sterowanie")
        hold off;

%         print("rozmyty_DMC.eps","-depsc","-r400")
    end
end