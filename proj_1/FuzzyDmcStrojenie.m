typ_funkcji = "gaus";
liczba_regulatorow = 5;

zaklocenia = 1; % 0 - brak zaklocen, 1 - zaklocenia
rysowanie = 0; % - brak rysowania wykresów w wywołaniu funkcji, 1 - rysowanie wykresów w wywołaniu funkcji
[E, h1, h2, h2_zad, F1, Fd]=FuzzyDmc([132, 1, 0.32], liczba_regulatorow, typ_funkcji, zaklocenia, rysowanie);
% do porównania
[Ed, h1d, h2d, h2_zadd, F1d, Fdd]=SuperDMC([132, 1, 0.32], zaklocenia, rysowanie);
disp(E)
disp(Ed)
figure(1)
hold on;
stairs(h2)
stairs(h2d)
stairs(h2_zad)
title("Wyjście")
legend("Wyjście rozmytego DMC","Wyjście konwencjonalnego DMC", "Wartość zadana")
xlabel("chwila k")
ylabel("wartość wyjścia")
hold off;
print("porownanieDMC_wyjscie_zak.eps","-depsc","-r400")

figure(2)
hold on;
stairs(F1)
stairs(F1d)
stairs(Fd)
legend("Sterowanie rozmytego DMC","Sterowanie konwencjonalnego DMC", "Zakłócenia", "Location", "SE")
xlabel("chwila k")
ylabel("wartość sterowania")
title("Sterowanie")
hold off;
print("porownanieDMC_sterowanie_zak.eps","-depsc","-r400")