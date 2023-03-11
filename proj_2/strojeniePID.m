clear all
close all
%Ograniczenia dolne parametrów
lb(1:3)=-200;
% %Ograniczenia górne parametrów
ub(1:3)=200;
%nastawy startowe
nastawy_startowe_PID = [-1, 10, 5];

wektor_PID = fmincon(@Pid,nastawy_startowe_PID, [], [], [], [], lb, ub);

% wektor_PID = [-0.606601881186912,3.53273449100721,3.71226835751801];
% wektor_PID = [-0.61, 3.52, 3.71]
[E, Y, yzad, U] = Pid(wektor_PID);

disp('PID')
disp(E)
figure(1)
stairs(Y)
hold on
stairs(yzad)
title(['Wyjście, błąd: ', num2str(E)])
legend("wyjście", "wartość zadana")
figure(2)
stairs(U)
title("Sterowanie")