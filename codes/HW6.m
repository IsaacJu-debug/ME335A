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

%% q2.0
close all 
k_array = [1,2,4];
i_array = [1,4,5,6];
nel_array = 2.^(i_array(:));

v_w_func = @(w, x) cos(w * x);
w_w_func = @(w, x) max(x.^w, 0);

u_func = @(x) v_w_func(30.0, x);
coord = [-1.0, 1.0];
plot_size = 50; % plot size per mesh
[int_func_array, plot_coord_array] = assemble(nel_array(4), k_array(2), coord, u_func, plot_size);
total_plot_size = plot_size * nel_array(4);
plot_func(plot_coord_array, int_func_array, u_func);

function plotConvergence(h_array, error_array)

    % Check inputs
    if length(h_array) ~= length(error_array)
        error('Input arrays must be of the same length.')
    end

    % Create figure
    figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7 5])

    % Create log-log plot
    loglog(h_array, error_array, 'ko-', 'LineWidth', 2, 'MarkerSize', 10)

    % Calculate and display the slope of the line
    P = polyfit(log(h_array(end-1:end)), log(error_array(end-1:end)), 1);
    disp(['The slope of the line is ', num2str(P(1))])

    % Add line showing slope
    hold on
    loglog(h_array, exp(P(2))*h_array.^P(1), 'r--', 'LineWidth', 2)

    % Add labels, title, and legend
    xlabel('Log(HMax)', 'FontSize', 14)
    ylabel('log(e_{a,2}(u_i - u_{i+1}))', 'FontSize', 14)
    title('Convergence Plot', 'FontSize', 16)
    legend('Data', ['Slope = ', num2str(P(1))], 'Location', 'best')
    grid on
    % Set axes properties
    set(gca, 'FontSize', 12, 'Box', 'on', 'XMinorTick', 'on', 'YMinorTick', 'on', 'XScale', 'log', 'YScale', 'log')

end

function plot_func(plot_coord_array, int_value_array, exact_func)
    % Input:
    % int_func - symbolic function with variable x
    % exact_func - anonymous function
    % coords - array of two endpoints

    y_values = exact_func(plot_coord_array);
    
    % Create a new figure
    figure;
    
    % Plot the symbolic function
    plot(plot_coord_array, int_value_array, 'LineWidth', 2, 'DisplayName', 'int\_func');
    hold on;
    
    % Plot the anonymous function
    plot(plot_coord_array, y_values, 'LineWidth', 2, 'DisplayName', 'exact\_func');
    
    % Add title and labels
    title('Comparison of int\_func and exact\_func');
    xlabel('x');
    ylabel('y');
    
    % Add legend
    legend('Location', 'best');
    
    hold off;
end


function [int_func_array, plot_coord_array] = assemble(nel, k, coord, u_func, plot_size)
    % nel: element number
    % k: polynomial order
    % coord: end points of the entire domain
    % u_func: ground truth
    %int_func_value = zeros(plot_size, 1);
    plot_coord_array = zeros(1, plot_size*nel );
    int_func_array = zeros(1, plot_size*nel );
    coord_array = linspace(coord(1), coord(2), nel + 1);
    for i = 1:nel
        % loop through element
        x_coord = coord_array(i: i+1);
        L = lagrange_interpolant(k, x_coord, u_func);
        plot_coord = linspace(x_coord(1), x_coord(2), plot_size);
        for j = 1:k+1
            % loop through local shape functions
            %disp(plot_coord(1,:));
            int_func_array(1, ( i-1)*plot_size + 1 : i*plot_size) = int_func_array(1, ( i-1)*plot_size + 1 : i*plot_size) +...
                                        double(subs(L(j), plot_coord(1,:)));
           
            plot_coord_array(1, ( i-1)*plot_size + 1 : i*plot_size) = plot_coord(1,:);
            %disp(double(subs(L(j), plot_coord(1,:))));
            %int_func = int_func + L(j);
        end 
    end
    
end

function L = lagrange_interpolant(k, x_coord, u_func)
    % Input: 
    % k - order of the Lagrange polynomial
    % x_coord - vector with start and end points of the interval [x_start, x_end]
    % u_func is the function to be approximated 

    % Create k+1 equally spaced nodes within the interval
    x_nodes = linspace(x_coord(1), x_coord(2), k+1);
    u_values = u_func(x_nodes);

    % Declare x as symbolic variable
    syms x;

    % Initialize L as a zero polynomial
    L = sym(zeros(1, k+1));

    % Compute the Lagrange polynomials
    for j = 1:(k+1)
        % Start with L_j(x) = 1
        L(j) = 1;

        % Multiply by (x - x_i) / (x_j - x_i) for all i ≠ j
        for i = 1:(k+1)
            if i ~= j
                L(j) = L(j) * (x - x_nodes(i)) / (x_nodes(j) - x_nodes(i));
            end
        end
        L(j) = L(j)*u_values(j);
    end

    %Modify the Lagrange polynomials so they are 0 outside of x_coord
    for j = 1:(k+1)
        L(j) = piecewise(x < x_coord(1), 0, x > x_coord(2), 0, L(j));
    end

end






