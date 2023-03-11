wektor = [130 1 0.32];
liczba_regulatorow = 5;
typ_funkcji = 'gaus';
zaklocenia = 1;
rysowanie = 0;


[E, h1, h2, h2_zad, F1, Fd] = SL(wektor, liczba_regulatorow, typ_funkcji, zaklocenia, rysowanie);
[En, h1n, h2n, h2_zadn, F1n, Fdn] = NPL(wektor, liczba_regulatorow, typ_funkcji, zaklocenia, rysowanie);
disp(E)
disp(En)

figure(1)
hold on
stairs(h2)
stairs(h2n)
stairs(h2_zad)
legend("Wyjście SL", "Wyjście NPL", "Wartość zadana")
print("porownanie_SL_NPL_wyjscie_zak.eps","-depsc","-r400")

figure(2)
hold on
stairs(F1)
stairs(F1n)
stairs(Fd)
legend("Sterowanie SL", "Sterowanie NPL", "Zakłócenia")
print("porownanie_SL_NPL_sterowanie_zak.eps","-depsc","-r400")