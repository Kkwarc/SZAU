function functions = MembershipFunction(liczba_regulatorow, typ_funkcji, y_min, y_max, rysowanie)
    width = y_max-y_min;
    functions = cell(1, liczba_regulatorow);
    for i = 1:1:liczba_regulatorow
        if typ_funkcji == "tr"
            functions{i} = ...
                @(x) trimf(x, [y_min+(i-1)*width/(liczba_regulatorow+1),...
                y_min+(i)*width/(liczba_regulatorow+1)...
                y_min+(i+1)*width/(liczba_regulatorow+1)]); 
        elseif typ_funkcji == "gaus"
            functions{i} =...
                @(x) gaussmf(x, [sqrt((width/(liczba_regulatorow+1))), ...
                y_min+i*width/(liczba_regulatorow+1)]);
        end
    end
    if rysowanie == 1
        x=y_min:0.01:y_max;
        figure(5)
        title("Funkcje przynależności")
        for i=1:1:liczba_regulatorow
            figure(5)
            plot(x, functions{i}(x))
            hold on
        end
        figure(5)
        print("funkcje_przynaleznosci_5_tr.eps","-depsc","-r400")
    end
end

