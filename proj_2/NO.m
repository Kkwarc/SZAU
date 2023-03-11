function [E, x] = NO(x, N, N_u, lambda, u, y, y_zad, k)
U = u;

U(k:k+N_u-1) = x;
U(k+N_u:k+N) = U(k+N_u-1);

dU = zeros(N_u, 1);
for i=0:N_u-1
    dU(k+i) = U(k+i)-U(k+i-1);
end

Y0(1) = siec(U(k-3), U(k-4), y(k), y(k-1));
Y0(2) = siec(U(k-2), U(k-3), Y0(1), y(k));
Y0(3) = siec(U(k-1), U(k-2), Y0(2), Y0(1));
for i=4:N
   Y0(i) = siec(U(k-4+i), U(k-5+i), Y0(i-1), Y0(i-2));
end

sum1 = 0;
sum2 = 0;
for i=1:N
    sum1 = sum1 + (y_zad(k)-Y0(i))^2;
end
for i=0:N_u-1
   sum2 = sum2 + dU(k+i)^2; 
end
E = sum1 +lambda*sum2;
end

