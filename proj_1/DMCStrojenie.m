opts = optimoptions('ga', 'MaxStallGenerations', 50, 'PopulationSize',20,"MaxGenerations",100);

down_lim = [1, 1, 0.01];
up_lim = [150, 100, 20];
calkowite_parametry = [1, 2];


% wektor = ga(@(x) SuperDMC(x, 0, 0), 3, [], [], [], [], down_lim, up_lim, [], calkowite_parametry, opts); % (x, 0) - 0 oznacza brak zaklocen
%  wektor = [132, 1, 0.32]; % parametry znalezione przez ga
wektor = [132, 1, 200]; % parametry znalezione poprawkami ręcznymi 

zaklocenia = 1; % 0 - brak zaklocen, 1 - zaklocenia
rysowanie = 0; % - brak rysowania wykresów w wywołaniu funkcji, 1 - rysowanie wykresów w wywołaniu funkcji
[E, h1, h2, h2_zad, F1, Fd]=SuperDMC(wektor, zaklocenia, rysowanie);

disp(E)
figure(1)
hold on
stairs(h2)
stairs(h2_zad)
legend("Wyjście", "Wartość zadana")
print("DMC_wyjscie_zak.eps","-depsc","-r400")
figure(2)
hold on
stairs(F1)
stairs(Fd)
legend("Sterowanie", "Zakłócenia")
print("DMC_sterowanie_zak.eps","-depsc","-r400")
