function u = losowe_sterowanie(liczbaKrokowDoZmiany, liczbaProbek, u_min, u_max)
    u(1:liczbaProbek) = 0;
    for k=2:(liczbaProbek/liczbaKrokowDoZmiany)-1
        u(k*liczbaKrokowDoZmiany:(k+1)*liczbaKrokowDoZmiany) = (u_max-u_min)*rand(1) + u_min;
    end
end

