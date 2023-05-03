function result = cal_bilinear_fun(u, w)
    % Define the symbolic variables
    syms x;

    % Calculate the symbolic derivatives of u and w
    u_prime = diff(u, x);
    w_prime = diff(w, x);

    % Define the integrand as a function of x
    integrand = (1 + x^2) * u_prime * w_prime + x * u_prime * w - x^2 * w * u;
    disp(integrand);
    % Calculate the integral of the integrand from 0 to 1
    result = int(integrand, x, 0, 1) - 6* subs(u, x, 1) * subs(w, x, 1);
    disp(result);
end

function result = cal_linear_fun(u, w)
    result = u(1)*w(1); 
end
