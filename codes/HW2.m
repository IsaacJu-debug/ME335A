%% on bases for vector spaces of functions
% Define the function N(x)
N = @(x, x0) max(1 - abs(x - x0), 0);

% Set the range for x
x = linspace(-3, 3, 1000);

% Calculate the function values for each x0
N_x0_minus1 = N(x, -1);
N_x0_0 = N(x, 0);
N_x0_1 = N(x, 1);

% Plot the functions
figure;
plot(x, N_x0_minus1, 'r-', 'LineWidth', 2); % red line
hold on;
plot(x, N_x0_0, 'g--', 'LineWidth', 2); % green dashed line
plot(x, N_x0_1, 'b-.', 'LineWidth', 2); % blue dash-dot line
hold off;

% Customize the plot
xlabel('x');
ylabel('N(x)');
title('Plot of N(x) with x0 = -1, x0 = 0, and x0 = 1');
legend('x0 = -1', 'x0 = 0', 'x0 = 1');
set(gca, 'Color', 'w');
grid on;
xlim([-3, 3]);
ylim([0, 1]);

%% 2
coeff_mat = [1/2, 1/2, 0, 1; 0, 1/2, 1/2, 1; 0, 1/4, 3/4, 1; 0, 0, 0, 0];
det(coeff_mat)


%% 3 
% Define the symbolic variable x
syms x;

% Define u_array and w_array
u_array = [x, x^2, x^3, 1];
w_array = [x, x^2, x^3, 1];

% Initialize k_mat as a 4x4 matrix of zeros
k_mat = sym(zeros(4, 4));
F = zeros(4,1);
% Calculate the entries of k_mat using bilinear_functional
for i = 1:length(u_array)
    for j = 1:length(w_array)
        if (i <= 3)
            k_mat(i, j) = cal_bilinear_fun(u_array(j), w_array(i));
        end
    end
end

k_mat(4, 4) = 1;
% Display k_mat
disp(k_mat);
disp(double(k_mat));
F(4) = 1;
u = double( (k_mat)^-1*F );
disp(u);

%%

syms u_x
u_x = u(1) * u_array(1) + u(2) * u_array(2) + u(3) * u_array(3) + u(4) * u_array(4);
u_x_prime = diff(u_x, x);
disp( double( 3 * subs(u_x, x, 1)));
disp( double(  subs(u_x_prime, x, 1)));
disp( double( 3 * subs(u_x, x, 1) - subs(u_x_prime, x, 1)));

%%
% Define the function u(x)
u = @(x) 1 - 1.5436*x + 0.0813*x.^2 + 0.1037*x.^3;

% Define the range of x values
x = linspace(0, 1, 100);

% Compute u(x) for each x value
y = u(x);

% Create a plot of u(x) against x
figure;
plot(x, y);
xlabel('x');
ylabel('u(x)');
title('Plot of u(x) = 1 - 1.5436x + 0.0813x^2 + 0.1037x^3');
grid on;

%% bilinear element 
% Define the range for x and y
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);

% Create a grid of points (X, Y) using the x and y vectors
[X, Y] = meshgrid(x, y);

% Compute the function values Z = f(X, Y) for the grid points
Z = 1/4 * (1 - X) .* (1 - Y);

% Plot the surface
surf(X, Y, Z);

% Add labels and a title
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('f(x, y) = 1/4 * (1 - x) * (1 - y)');

