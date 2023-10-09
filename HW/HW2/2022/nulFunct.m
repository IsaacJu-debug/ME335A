function [outcome] = nulFunct(v, x)
    outcome = int(- v*(x^2), x, 0, 1) - 6;
    return 
end