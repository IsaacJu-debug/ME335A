function [outcome] = iLegendre(n, x)
    res = 0;
    coeff = ((-1)^(n+1))*sqrt(2*n + 1);
    for k = 0:n
        res = res + (nchoosek(n,k) * nchoosek(n+k,k) * (-x)^(k+1))/(k+1);
    end
    outcome = res*coeff;
    return
end