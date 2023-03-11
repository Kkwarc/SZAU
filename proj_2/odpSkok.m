function s = odpSkok()
    u(1:10) = 0;
    u(11:100) = 1;
    
    alpha1 = -1.599028;
    alpha2 = 0.632337;
    betha1 = 0.010754;
    betha2 = 0.009231;
    
    x1(1:100) = 0;
    x2(1:100) = 0;
    y(1:100) = 0;
    for k=10:100
        g_1 = g1(u(k-4));
        x1(k) = -alpha1*x1(k-1)+x2(k-1)+betha1*g_1;
        x2(k) = -alpha2*x1(k-1)+betha2*g_1;
        y(k) = g2(x1(k));
    end
    s = y(11:end);
    stairs(s)
end

