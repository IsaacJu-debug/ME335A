function [outcome] = alFunct(u, v, x)
    outcome = int((1 + x^2)*diff(u)*diff(v) + x*v*diff(u) - u*v*(x^2), x, 0, 1) - 6;
    return 
end