%% q1.1.1 plot functions
close all
% Define x in the interval (0, pi)
x = linspace(0,pi,1000);

% Define n values
n_values = 2.^(1:5);

% Create and plot f_n(x)
figure
hold on
for n = n_values
    fn = n * sin(x./n);
    plot(x,fn,'DisplayName',sprintf('n = %d',n))
end
title('Plot of f_n(x) = n * sin(x / n)')
xlabel('x')
ylabel('f_n(x)')
legend('show','Location','best')
hold off

% Create and plot g_n(x)
figure
hold on
for n = n_values
    gn = 1./n * sin(n.*x);
    plot(x,gn,'DisplayName',sprintf('n = %d',n))
end
title('Plot of g_n(x) = 1/n * sin(nx)')
xlabel('x')
ylabel('g_n(x)')
legend('show','Location','best')
hold off

% Create and plot h_n(x)
figure
hold on
for n = n_values
    hn = 1./(1 + n.*x);
    plot(x,hn,'DisplayName',sprintf('n = %d',n))
end
title('Plot of h_n(x) = 1 / (1 + nx)')
xlabel('x')
ylabel('h_n(x)')
legend('show','Location','best')
hold off

%% q1.1.2 plot derivatives 

% Define x in the interval (0, pi)
x = linspace(0,pi,1000);
dx = x(2) - x(1);  % compute delta x for derivative

% Define n values
n_values = 2.^(1:5);

% Create and plot derivatives of f_n(x)
figure
hold on
for n = n_values
    fn = n * sin(x./n);
    dfn = diff(fn)./dx;
    plot(x(1:end-1),dfn,'DisplayName',sprintf('n = %d',n))
end
title('Plot of df_n/dx')
xlabel('x')
ylabel('df_n/dx')
legend('show','Location','best')
hold off

% Create and plot derivatives of g_n(x)
figure
hold on
for n = n_values
    gn = 1./n * sin(n.*x);
    dgn = diff(gn)./dx;
    plot(x(1:end-1),dgn,'DisplayName',sprintf('n = %d',n))
end
title('Plot of dg_n/dx')
xlabel('x')
ylabel('dg_n/dx')
legend('show','Location','best')
hold off

% Create and plot derivatives of h_n(x)
figure
hold on
for n = n_values
    hn = 1./(1 + n.*x);
    dhn = diff(hn)./dx;
    plot(x(1:end-1),dhn,'DisplayName',sprintf('n = %d',n))
end
title('Plot of dh_n/dx')
xlabel('x')
ylabel('dh_n/dx')
legend('show','Location','best')
hold off

%% q1.3 
syms x n real;  % declare x and n as symbolic variables

f = n * sin(x/n);  % define the function
f_squared = f^2;  % square the function

% compute the integral
integral_f_squared = int(f_squared, x, 0, pi);

% display the result
disp(integral_f_squared);

f_derivative = diff(f, x);  % compute the derivative of f
disp(f_derivative);
f_derivative_squared = f_derivative^2;  % square the derivative

% compute the integral of the squared derivative
integral_f_derivative_squared = int(f_derivative_squared, x, 0, pi);

% display the result
disp(integral_f_derivative_squared);

%% q1.3.2
close all
% define f_l2 and f_h1
f_l2 = integral_f_squared^0.5;
f_h1 = (integral_f_derivative_squared + integral_f_squared)^0.5;

% create an array of n values
n_values = 1:0.1:10;

% convert symbolic functions to function handles for numerical evaluation
f_l2_handle = matlabFunction(f_l2);
f_h1_handle = matlabFunction(f_h1);

% compute f_l2 and f_h1 values for given n values
f_l2_values = f_l2_handle(n_values);
f_h1_values = f_h1_handle(n_values);

% create and configure the plot
figure
plot(n_values, f_l2_values, 'DisplayName', '||f||_{0,2}')
hold on
plot(n_values, f_h1_values, 'DisplayName', '||f||_{1,2}')
xlabel('n')
ylabel('Function Value')
title('Plot of ||f||_{0,2} and ||f||_{1,2} vs n')
legend('show','Location','best')
hold off
%% q1.4.1
close all
% declare x and n as symbolic variables
syms x n real;

% define the functions
f_n = n * sin(x/n); 
f_inf = x; 

% calculate the difference
f_diff = f_n - f_inf;

% calculate the square of the difference
f_diff_squared = f_diff^2;

% calculate the square of the derivative of the difference
f_diff_derivative_squared = (diff(f_diff, x))^2;

% compute the integrals of the squares and their square roots
l2_norm = sqrt(int(f_diff_squared, x, 0, pi));
h1_norm = sqrt(int(f_diff_squared + f_diff_derivative_squared, x, 0, pi));

% convert symbolic expressions to function handles for numerical evaluation
l2_norm_handle = matlabFunction(l2_norm, 'Vars', {n});
h1_norm_handle = matlabFunction(h1_norm, 'Vars', {n});
f_diff_handle = matlabFunction(f_diff, 'Vars', {n, x});

% create an array of n values
n_values = linspace(1, 50, 100);
x = linspace(0, pi, 50);
% compute l2_norm and h1_norm values for given n values
l2_norm_values = arrayfun(l2_norm_handle, n_values);
h1_norm_values = arrayfun(h1_norm_handle, n_values);

% compute inf_norm values for each n using numerical optimization
inf_norm_values = zeros(size(n_values));

for i = 1:numel(n_values)
    inf_norm_values(i) = max(abs(f_diff_handle(n_values(i), x)));
end

% create and configure the plot
% figure
% hold on
% plot(n_values, l2_norm_values, 'DisplayName', 'L2 Norm')
% plot(n_values, inf_norm_values, 'DisplayName', 'Infinity Norm')
% plot(n_values, h1_norm_values, 'DisplayName', 'H1 Norm')
% xlabel('n')
% ylabel('Norm')
% title('Plot of L2, Infinity, and H1 Norms vs n')
% legend('show','Location','best')
% hold off
%%
% create and configure the plot
figure
loglog(n_values, l2_norm_values, 'LineWidth', 2, 'DisplayName', 'L2 Norm')
hold on
loglog(n_values, inf_norm_values, 'LineWidth', 2, 'DisplayName', 'Infinity Norm')
loglog(n_values, h1_norm_values, 'LineWidth', 2, 'DisplayName', 'H1 Norm')
xlabel('n')
ylabel('Norm')
title('Log-Log Plot of L2, Infinity, and H1 Norms vs n')
legend('show','Location','best')
grid on
hold off

%%
close all
% declare x and n as symbolic variables
syms x n real;

% define the functions
f_n = 1/n * sin(n*x); 
f_inf = 0; 

% calculate the difference
f_diff = f_n - f_inf;

% calculate the square of the difference
f_diff_squared = f_diff^2;

% calculate the square of the derivative of the difference
f_diff_derivative_squared = (diff(f_diff, x))^2;

% compute the integrals of the squares and their square roots
l2_norm = sqrt(int(f_diff_squared, x, 0, pi));
h1_norm = sqrt(int(f_diff_squared + f_diff_derivative_squared, x, 0, pi));

% convert symbolic expressions to function handles for numerical evaluation
l2_norm_handle = matlabFunction(l2_norm, 'Vars', {n});
h1_norm_handle = matlabFunction(h1_norm, 'Vars', {n});
f_diff_handle = matlabFunction(f_diff, 'Vars', {n, x});

% create an array of n values
n_values = linspace(1, 50, 100);
x = linspace(0, pi, 50);
% compute l2_norm and h1_norm values for given n values
l2_norm_values = arrayfun(l2_norm_handle, n_values);
h1_norm_values = arrayfun(h1_norm_handle, n_values);

% compute inf_norm values for each n using numerical optimization
inf_norm_values = zeros(size(n_values));

for i = 1:numel(n_values)
    inf_norm_values(i) = max(abs(f_diff_handle(n_values(i), x)));
end

% create and configure the plot
% figure
% hold on
% plot(n_values, l2_norm_values, 'DisplayName', 'L2 Norm')
% plot(n_values, inf_norm_values, 'DisplayName', 'Infinity Norm')
% plot(n_values, h1_norm_values, 'DisplayName', 'H1 Norm')
% xlabel('n')
% ylabel('Norm')
% title('Plot of L2, Infinity, and H1 Norms vs n')
% legend('show','Location','best')
% hold off
figure
loglog(n_values, l2_norm_values, 'LineWidth', 2, 'DisplayName', 'L2 Norm')
hold on
loglog(n_values, inf_norm_values, 'LineWidth', 2, 'DisplayName', 'Infinity Norm')
loglog(n_values, h1_norm_values, 'LineWidth', 2, 'DisplayName', 'H1 Norm')
xlabel('n')
ylabel('Norm')
title('Log-Log Plot of L2, Infinity, and H1 Norms vs n')
legend('show','Location','best')
grid on
hold off

%%
close all
% declare x and n as symbolic variables
syms x n real;

% define the functions
f_n = 1/(1 + n *x); 
f_inf = 0; 

% calculate the difference
f_diff = f_n - f_inf;

% calculate the square of the difference
f_diff_squared = f_diff^2;

% calculate the square of the derivative of the difference
f_diff_derivative_squared = (diff(f_diff, x))^2;

% compute the integrals of the squares and their square roots
%l2_norm = sqrt(int(f_diff_squared, x, 0.0, pi));
l2_norm = (pi / (n * pi + 1))^0.5;
%h1_norm = sqrt(int(f_diff_squared + f_diff_derivative_squared, x, 0.0, pi));
h1_norm = ((n^2/3 + 1)/n - ((pi^2 + 1/3)*n^2 + 2*pi*n + 1)/(n*(n*pi + 1)^3)) ^0.5; 

% convert symbolic expressions to function handles for numerical evaluation
l2_norm_handle = matlabFunction(l2_norm, 'Vars', {n});
h1_norm_handle = matlabFunction(h1_norm, 'Vars', {n});
f_diff_handle = matlabFunction(f_diff, 'Vars', {n, x});

% create an array of n values
n_values = linspace(1, 50, 100);
x = linspace(0, pi, 50);
% compute l2_norm and h1_norm values for given n values
l2_norm_values = arrayfun(l2_norm_handle, n_values);
h1_norm_values = arrayfun(h1_norm_handle, n_values);

% compute inf_norm values for each n using numerical optimization
inf_norm_values = zeros(size(n_values));

for i = 1:numel(n_values)
    inf_norm_values(i) = max(abs(f_diff_handle(n_values(i), x)));
end

% create and configure the plot
% figure
% hold on
% plot(n_values, l2_norm_values, 'DisplayName', 'L2 Norm')
% plot(n_values, inf_norm_values, 'DisplayName', 'Infinity Norm')
% plot(n_values, h1_norm_values, 'DisplayName', 'H1 Norm')
% xlabel('n')
% ylabel('Norm')
% title('Plot of L2, Infinity, and H1 Norms vs n')
% legend('show','Location','best')
% hold off
figure
loglog(n_values, l2_norm_values, 'LineWidth', 2, 'DisplayName', 'L2 Norm')
hold on
loglog(n_values, inf_norm_values, 'LineWidth', 2, 'DisplayName', 'Infinity Norm')
loglog(n_values, h1_norm_values, 'LineWidth', 2, 'DisplayName', 'H1 Norm')
xlabel('n')
ylabel('Norm')
title('Log-Log Plot of L2, Infinity, and H1 Norms vs n')
legend('show','Location','best')
grid on
hold off

%%
