function [E, h1, h2, h2_zad, F1, Fd]=SL(wektor, liczba_regulatorow, typ_funkcji, zaklocenia, rysowanie)
% wektor = [130, 1 , 0.32]
% liczba_regulatorow =5;
% typ_funkcji = 'gaus';
% zaklocenia = 0;
% rysowanie =1;
E=0;

T=5;
tau = round(160/T);

y_min = 5;
y_max = 35;

odp_skok = StepResponsesFuzzy(liczba_regulatorow, [y_min, y_max], 0);

functions = MembershipFunction(liczba_regulatorow, typ_funkcji, y_min, y_max, 0);

% 
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
poczatek = start; %chwila k w której zmienia sie wartość zadana zalecana > 80

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
        sum_mi = sum_mi + max(functions{i}(h2(k)), 0.01);
    end
    
    s_average = zeros(1, D);
    for j=1:D
        for i=1:liczba_regulatorow
            s_average(j) = s_average(j) + max(functions{i}(h2(k)), 0.01)*odp_skok{i}(j)/sum_mi;
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

    M_p = zeros(N, D-1);
    for column=1:(D-1)
        for row=1:N
            if row + column > D
                if column>D
                    M_p(row, column) = 0;
                else
                    M_p(row, column) = s_average(D) - s_average(column);
                end
            else
                M_p(row, column) = s_average(row + column) - s_average(column);
            end
        end
    end


    %Obliczenie DU_p
    for d=1:(D-1)
        DU_p(d) = F1(k-d) - F1(k-d-1);
    end

    %Obliczenie sterowania

    %Pomiar wyjścia
    Y = ones(N, 1) * h2(k);

    Y_zad = ones(N, 1) * h2_zad(k);

    Y0 =  M_p * DU_p + Y;

%     cel = @(du) (Y_zad-Y0 - M*du)'*(Y_zad-Y0 - M*du)+lambda*du'*du;

    
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
        figure(2)
        hold on
        stairs(F1)
        title("Sterowanie")
    end
    
end

function [dh1, dh2] = func(h1, h2, F1, Fd)
    global C1 alpha1 C2 alpha2 
    dh1 = (1/(3*C1*h1^2)) * (F1 + Fd - alpha1*sqrt(h1));
    dh2 =  (1/(3*C2*h2^2)) * (alpha1*sqrt(h1) - alpha2*sqrt(h2));
end