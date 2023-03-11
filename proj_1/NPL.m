function [E, h1, h2, h2_zad, F1, Fd]=NPL(wektor, liczba_regulatorow, typ_funkcji, zaklocenia, rysowanie)
global C1 alpha1 C2 alpha2 Fd

E=0;

T=5;
tau = round(160/T);

y_min = 5;
y_max = 35;

odp_skok = StepResponsesFuzzy(liczba_regulatorow, [y_min, y_max], 0);

functions = MembershipFunction(liczba_regulatorow, typ_funkcji, y_min, y_max, 0);

x=y_min:0.01:y_max;
figure(3)
for i=1:1:liczba_regulatorow
    title("Funkcje przynależności")
    plot(x, functions{i}(x))
    hold on
end

D = 1000;

start= D+1;

endt=5000+D; %koniec symulacji
poczatek = start; %chwila k w której zmienia sie wartość zadana

%Definicja horyzontów i parametrów
N = wektor(1);

N_u = wektor(2);

lambda = wektor(3);

%deklaracja potrzebnych wektorów
h1(1:1:start)=18.9225;
h1(start+1:1:endt)=0;
h2(1:1:start)=15.6548;
h2(start+1:1:endt)=0;

C1=0.35;
C2 = 0.3;
alpha1 = 20;
alpha2 = 22;
Fd = 15;

u_min = sqrt(y_min)*alpha2-15;
u_max = sqrt(y_max)*alpha2-15;

F1(1:1:start+1) = 72;
F1(start+2:1:endt) = 0;

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

options = optimoptions('quadprog','Display','off');

for k=start:endt
    disp("Chwila: "+(k-start)+"/"+(endt-start))
   
    %symulacja obiektu
    [h1_prov, h2_prov] = func(h1(k-1), h2(k-1), F1(k-1-tau), Fd(k));
    [h1_mid, h2_mid] = func(h1(k-1)+h1_prov*T/2, h2(k-1)+h2_prov*T/2, F1(k-1-tau), Fd(k));
    h1(k) = h1(k-1) + T*h1_mid ;
    h2(k) = h2(k-1) + T*h2_mid;

    %Obliczenie części macierzy DMC

    sum_mi = 0;
    for i=1:liczba_regulatorow
        sum_mi = sum_mi + functions{i}(h2(k));
    end
    
    s_average = zeros(1, D);
    for j=1:D
        for i=1:liczba_regulatorow
            s_average(j) = s_average(j) + functions{i}(h2(k))*odp_skok{i}(j)/sum_mi;
        end
    end

    M = zeros(N, N_u);
    for column=1:N_u
        for row=1:N
           if (row>=column)
             if(row-column+1<=D)
                M(row,column)=s_average(row-column+1);
             else
               M(row,column)=s_average(D);
             end
          end
        end
    end
    %Obliczenie DU_p
    for d=1:(D-1)
        DU_p(d) = F1(k-d) - F1(k-d-1);
    end

    %Obliczenie sterowania

    %Pomiar wyjścia

    Y_zad = ones(N, 1) * h2_zad(k);
    d = h2(k) - s_average(1:D-1)*DU_p - s_average(D)*F1(k-D);

    Y0 =  zeros(N,1);
    h01 = zeros(N,1);

    s_av = s_average(2:D); 

    
    for i=1:N
        if i>1
            [h01_prov, Y0_prov] = func(h01(i-1), Y0(i-1), F1(min(k+i-1-tau,k-1)), Fd(k));
            [h01_mid, Y0_mid] = func(h01(i-1)+h01_prov*T/2, Y0(i-1)+Y0_prov*T/2, F1(min(k+i-1-tau,k-1)),Fd(k));
            h01(i) = h01(i-1) + T*h01_mid ;
            Y0(i) = Y0(i-1) + T*Y0_mid;
        else
            [h01_prov, Y0_prov] = func(h1(k), h2(k), F1(k-tau), Fd(k));
            [h01_mid, Y0_mid] = func(h1(k)+h01_prov*T/2, h2(k)+Y0_prov*T/2, F1(k-tau),Fd(k));
            h01(i) = h1(k) + T*h01_mid ;
            Y0(i) = h2(k) + T*Y0_mid;
        end
        Y0(i) = max(y_min + 0.001, Y0(i));
        Y0(i) = min(y_max - 0.001, Y0(i));
    end

    
    A = [-1*ones(1, N_u); 1*ones(1,N_u)]; 
    B = [F1(k-1)-u_min-0.001;u_max-0.001-F1(k-1)];

    DU = quadprog(2*(M')*M+2*lambda*eye(N_u,N_u), -2*(Y_zad-Y0)'*M,A,B,[],[],[],[],[], options);

    F1(k) = F1(k-1) + DU(1);

    E=E+(h2_zad(k)-h2(k))^2;

end
    if rysowanie == 1
        figure(1)
        hold on
        stairs(h2)
        stairs(h2_zad)
        title("Wyjście")
%         print("wyjscie_4.eps","-dpng", "-r400")
        figure(2)
        hold on
        stairs(F1)
        stairs(Fd)
        legend("sterowanie", "zakłócenia")

%         print("wejscie_4.eps","-dpng", "-r400")
    end
end

function [dh1, dh2] = func(h1, h2, F1, Fd)
    global C1 alpha1 C2 alpha2 
    dh1 = (1/(3*C1*h1^2)) * (F1 + Fd - alpha1*sqrt(h1));
    dh2 =  (1/(3*C2*h2^2)) * (alpha1*sqrt(h1) - alpha2*sqrt(h2));
end
