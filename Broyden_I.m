function [x, k] = Broyden_I(x0, A, f1, f2, k)

    for i = 1 : k

       x1 = x0 - inv(A) * F(f1, f2, x0);

       r = x1 - x0;

       tri = F(f1, f2,x1) - F(f1, f2, x0);

       A = A + ((tri - A * r) * r') / (r' * r);

       x0 = x1;

    end

    x = x1;

end