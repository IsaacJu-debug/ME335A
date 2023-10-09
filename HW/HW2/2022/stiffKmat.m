function [K] = stiffKmat(x)
    K = zeros(4, 4);
    %K(1,1) = alFunct(x,x,x);
    %K(1,2) = alFunct(x^2, x, x);
    %K(1,3) = alFunct(x^3,x,x);
    %K(1,4) = alFunct(1,x,x);
    K
    K = [alFunct(x,x,x) alFunct(x^2,x,x) alFunct(x^3,x,x) nulFunct(x,x);
        alFunct(x, x^2, x) alFunct(x^2, x^2, x) alFunct(x^3, x^2, x) nulFunct(x^2, x);
        alFunct(x, x^3, x) alFunct(x^2, x^3, x) alFunct(x^3, x^3, x) nulFunct(x^3, x);
        0 0 0 1];
    return 
end