function [outcome] = integral_F(n, x)
    outcome = int((cos(4*pi*x) * iLegendre(n, x)), 0, 1);
    return
end