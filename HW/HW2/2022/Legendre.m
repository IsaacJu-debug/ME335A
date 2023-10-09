function [outcome] = Legendre(n, x)
    res = 0;
    coeff = ((-1)^n)*sqrt(2*n + 1);
    for k = 0:n
        res = res + (nchoosek(n,k) * nchoosek(n+k,k) * (-x)^k);
    end
    outcome = res*coeff;
    return
end