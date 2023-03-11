% function [E, x1, x2, y_zad, u]=NPL(wektor,rysowanie)
clear all
wektor = [50, 5, 5];
rysowanie = 1;
disp(wektor)

E=0;

% odp_skok = StepResponsesFuzzy(liczba_regulatorow, [y_min, y_max], 0);

D = 100;

start= D+1;

endt=600+D; %koniec symulacji
poczatek = start; %chwila k w której zmienia sie wartość zadana

%Definicja horyzontów i parametrów
N = wektor(1);

N_u = wektor(2);

lambda = wektor(3);

%deklaracja potrzebnych wektorów
x1(1:endt)=0;
x2(1:endt)=0;
u(1:endt)=0;
y(1:endt)=0;

alpha1 = -1.599028;
alpha2 = 0.632337;
betha1 = 0.010754;
betha2 = 0.009231;

%skok wartości zadanej
y_zad(1:poczatek)=0; 
y_zad(poczatek+1:150+D)=2;
y_zad(151+D:250+D)=-1;
y_zad(251+D:350+D)=3;
y_zad(351+D:endt)= 0;

reg = 'PID'; %'NPL', 'GPC', 'PID', 'NO'


if strcmp(reg, 'PID')
    K=-0.61; % człon proporcjonalny
    Ti=3.52; % człon całkujący
    Td=3.71; % człon różniczkujący 
    
    T = 1;
    
    r2 = (K*Td)/T;
    r1 = K*(T/(2*Ti) - 2*(Td/T) - 1);
    r0 = K*(1+ T/(2*Ti) + Td/T);
    
    e(1:endt)=0;
end

if strcmp(reg, 'GPC')
    S = zeros(1, N);
    naj_kward = [-0.0624180043648385;0.0512091664509879;-1.93950946976507;0.942753195786931];
    b4 = naj_kward(1);
    b5 = naj_kward(2);
    a1 = naj_kward(3);
    a2 = naj_kward(4);
    b = [0, 0, 0, b4, b5];
    a = [a1, a2];
    nb = 5;
    na = 2;
    for p=1:N
        s = 0;
        for i=1:min(p, nb)
            s = s + b(i);
        end
        for i=1:min(p-1, na)
            s = s - a(i)*S(p-i);
        end
        S(p) = s;
    end
end
    
for k=start:endt
%     disp("Chwila: "+(k-start)+"/"+(endt-start))
   
    %symulacja obiektu
    g_1 = g1(u(k-4));
    x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
    x2(k) = -alpha2*x1(k-1)+betha2*g_1;
    y(k) = g2(x1(k));
    
    if strcmp(reg, 'NO')
        u0 = u(k-1)*ones(1, N_u);
        opts = optimset('Display', 'off');
        disp(endt-k)
        U = fmincon(@(x) NO(x, N, N_u, lambda, u, y, y_zad, k), u0, [], [], [], [], -1*ones(1, N_u), 1*ones(1, N_u), [], opts);
        u(k) = U(1);
    elseif strcmp(reg, 'PID')
        e(k)=y_zad(k)-y(k);
        u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    else
        Yzad(1:N,1) = y_zad(k);

        if strcmp(reg, 'NPL')
            theta = 1*10^-5;

            f0 = siec(u(k-4), u(k-5), y(k-1), y(k-2));
            b4 = (siec(u(k-4)+theta, u(k-5), y(k-1), y(k-2))-f0)/theta;
            b5 = (siec(u(k-4), u(k-5)+theta, y(k-1), y(k-2))-f0)/theta;
            a1 = -(siec(u(k-4), u(k-5), y(k-1)+theta, y(k-2))-f0)/theta;
            a2 = -(siec(u(k-4), u(k-5), y(k-1), y(k-2)+theta)-f0)/theta;
            b = [0, 0, 0, b4, b5];
            a = [a1, a2, 0, 0];

            S = zeros(1, N);
            nb = 5;
            na = 2;
            for p=1:N
                s = 0;
                for i=1:min(p, nb)
                    s = s + b(i);
                end
                for i=1:min(p-1, na)
                    s = s - a(i)*S(p-i);
                end
                S(p) = s;
            end  
        end

        M = zeros(N, N_u);
        for row=1:N
            for column=1:N_u
                if row >= column
                    M(row, column) = S(row-column+1);
                end
            end
        end

        K = (M'*M+lambda*eye(N_u, N_u))^(-1)*M';

        Y0 = zeros(N, 1);
        if strcmp(reg, 'NPL')
            d = y(k) - f0;
            Y0(1) = siec(u(k-3), u(k-4), y(k), y(k-1))+d;
            Y0(2) = siec(u(k-2), u(k-3), Y0(1), y(k))+d;
            Y0(3) = siec(u(k-1), u(k-2), Y0(2), Y0(1))+d;
            for i=4:N
                Y0(i) = siec(u(k-1), u(k-1), Y0(i-1), Y0(i-2))+d;
            end
        end
        if strcmp(reg, 'GPC')
            f0 = b4*u(k-4)+b5*u(k-5)-a1*y(k-1)-a2*y(k-2);
            d = y(k) - f0;
            Y0(1) = b4*u(k-3)+b5*u(k-4)-a1*y(k)-a2*y(k-1)+d;
            Y0(2) = b4*u(k-2)+b5*u(k-3)-a1*Y0(1)-a2*y(k)+d;
            Y0(3) = b4*u(k-1)+b5*u(k-2)-a1*Y0(2)-a2*Y0(1)+d;
            for i=4:N
                Y0(i) = b4*u(k-1)+b5*u(k-1)-a1*Y0(i-1)-a2*Y0(i-2)+d;
            end
        end

        DU = K*(Yzad-Y0);

        u(k) = u(k-1) + DU(1);
    end
    if u(k) > 1
        u(k) = 1;
    end
    if u(k) < -1
        u(k) = -1;
    end
    E=E+(y_zad(k)-y(k))^2;

end
disp(reg)
disp(E)
    if rysowanie == 1
        figure(1)
        subplot(2,1,1)
        hold on
        stairs(y)
        stairs(y_zad)
        legend('Wyjście', 'Wartość zadana')
        title(["Wyjście"+", Błąd: "+ num2str(E)])
        subplot(2,1,2)
        stairs(u)
        title("Sterowanie")
        print("PID.eps","-depsc","-r400")
    end
% end

