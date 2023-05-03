%% on Convergence
% Note the following symbolic calculation is not accurate
syms x
syms u(x)
syms f(x)
syms F(x)
syms u_prime
u(x) = -(sin(2*pi*x))^2 / (8*pi^2);
f(x) = cos(4*pi* x);
u_prime = diff(u(x), x);

F(x) = diff( diff(u(x), x), x) + f(x);
disp(F(x))

u_prime
subs(u_prime, x, 1)

%% Finite element approximation

% Define domain and intervals
x = linspace(0, 1, 1000);
interval = 0.25;

% Initialize the plot
figure;
hold on;

% Create a cell array to store the legend strings
legend_strings = cell(1, 5);

% Define the hat functions
for i = 0:4
    % Define the left, middle and right points for each hat function
    left = max(0, i * interval - interval);
    middle = i * interval;
    right = min(1, i * interval + interval);
    
    % Define the piece-wise function
    y = zeros(size(x));
    y(x >= left & x < middle) = (x(x >= left & x < middle) - left) / (middle - left);
    y(x >= middle & x <= right) = (right - x(x >= middle & x <= right)) / (right - middle);
    
    % Plot the hat function
    plot(x, y);
    
    % Add the legend string for the current hat function
    legend_strings{i+1} = sprintf('N_%d', i + 1);
end

% Customize the plot
title('Piece-wise Hat Functions');
xlabel('x');
ylabel('y');
xlim([0, 1]);
ylim([0, 1.0]);
grid on;
legend(legend_strings);
hold off;

%% Stiffness matrix
syms x

k21_int = -1/0.25 * 1/0.25 + 2 * x/0.25 * (0.25 - x)/0.25;
k21 = int(k21_int, 0, 0.25);
k21 

k22_int_part1 = 1/0.25 * 1/0.25 + 2 * ((x)/0.25)^2;
k22_int_part2 = 1/0.25 * 1/0.25 + 2 *((0.5 - x)/0.25)^2;

k22= int(k22_int_part1, 0, 0.25) + int(k22_int_part2, 0.25, 0.5);
k22

k23_int = -1/0.25 * 1/0.25 + 2 * (x - 0.25)/0.25 * (0.5 - x)/0.25;
k23 = int(k23_int, 0.25, 0.5);
k23

k33_int_part1 = 1/0.25 * 1/0.25 + 2 * ((x -0.25)/0.25)^2;
k33_int_part2 = 1/0.25 * 1/0.25 + 2 *((0.75 - x)/0.25)^2;
k33= int(k33_int_part1, 0.25, 0.5) + int(k33_int_part2, 0.5, 0.75);
k33

k34_int = -1/0.25 * 1/0.25 + 2 * (x - 0.5)/0.25 * (0.75 - x)/0.25;
k34 = int(k34_int, 0.5, 0.75);
k34


k55_int = 1/0.25 * 1/0.25 +  2 * ((x -0.75)/0.25)^2;
k55 = int(k55_int, 0.75, 1.0) + 1;
k55

%% Load vector
syms x

N2_int_part1 = x/0.25 * x^2;
N2_int_part2 = ((0.5 - x)/0.25) * x^2;
N2_load = int(N2_int_part1, 0, 0.25) + int(N2_int_part2, 0.25, 0.5);
N2_load

N3_int_part1 = (x-0.25)/0.25 * x^2;
N3_int_part2 = ((0.75 - x)/0.25) * x^2;
N3_load = int(N3_int_part1, 0.25, 0.5) + int(N3_int_part2, 0.5, 0.75);
N3_load

N4_int_part1 = (x-0.5)/0.25 * x^2;
N4_int_part2 = ((1.0 - x)/0.25) * x^2;
N4_load = int(N4_int_part1, 0.5, 0.75) + int(N4_int_part2, 0.75, 1.0);
N4_load

N5_int = (x - 0.75) / 0.25 * x^2;
N5_load = int(N5_int, 0.75, 1.0) + 1;
N5_load

%% solve the linear system
K = [  1,      0,      0,      0,      0;
     -47/12,  25/3,  -47/12,   0,      0;
       0,    -47/12,  25/3,  -47/12,   0;
       0,      0,    -47/12,  25/3,  -47/12;
       0,      0,      0,    -47/12,  31/6];

F = [0, 7/384, 25/384, 55/384, 283/256];
U = K\F';
U
