
%% Q1-1 Lagrange p-k element

% Define the range and the number of points for plotting
x = linspace(0, 1, 1000);

% Calculate the Lagrange P-3 basis functions within the range [0, 1]
[L0, L1, L2, L3] = lagrange_p3_basis(x, 0, 1);

% Calculate N1 as the sum of L0, L1, L2, and L3
N1 = L0 + L1 + L2 + L3;

% Plot the Lagrange P-3 basis functions
figure;
plot(x, L0, 'r', 'LineWidth', 2);
hold on;
plot(x, L1, 'g', 'LineWidth', 2);
plot(x, L2, 'b', 'LineWidth', 2);
plot(x, L3, 'c', 'LineWidth', 2);
plot(x, N1, '--k', 'LineWidth', 2);

% Configure the plot
title('Lagrange P-3 Shape Functions in element 1');
xlabel('x');
ylabel('Shape Functions');
legend({'L0', 'L1', 'L2', 'L3', 'N1'}, 'Location', 'Best');
grid on;

%% Q1-1-1
% Define the range and the number of points for plotting
x = linspace(0, 3, 1000);

% Define N1 for each subplot
N1a = (x >= 0 & x <= 1);
N1b = (x >= 1 & x <= 2);
N1c = (x >= 2 & x <= 3);

% Create the subplots
figure;

% First subplot
subplot(1, 3, 1);
plot(x, N1a, 'r', 'LineWidth', 2);
ylim([-0.1, 1.1]);
title('N1(x)');
xlabel('x');
ylabel('N1(x)');
grid on;

% Second subplot
subplot(1, 3, 2);
plot(x, N1b, 'g', 'LineWidth', 2);
ylim([-0.1, 1.1]);
title('N2(x)');
xlabel('x');
ylabel('N2(x)');
grid on;

% Third subplot
subplot(1, 3, 3);
plot(x, N1c, 'b', 'LineWidth', 2);
ylim([-0.1, 1.1]);
title('N3(x)');
xlabel('x');
ylabel('N3(x)');
grid on;

%% Q1-1-2
% Define the range and the number of points for plotting
x1 = linspace(0, 1, 100);
x2 = linspace(1, 2, 100);
x3 = linspace(2, 3, 100);

% Calculate the Lagrange P-3 basis functions within the ranges [0, 1], [1, 2], and [2, 3]
[L10, L11, L12, L13] = lagrange_p3_basis(x1, 0, 1);
[L20, L21, L22, L23] = lagrange_p3_basis(x2, 1, 2);
[L30, L31, L32, L33] = lagrange_p3_basis(x3, 2, 3);

C1 = zeros(100);
C2 = zeros(100);
C3 = zeros(100);

% Combine x1, x2, and x3 for plotting
x = [x1, x2, x3];

% Create the subplots
figure;

% First subplot
subplot(2, 3, 1);
plot(x1, L10, 'r', x2, L20, 'g', x3, L30, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N1(x)');
xlabel('x');
ylabel('N1(x)');
grid on;

% Second subplot
subplot(2, 3, 2);
plot(x1, L11, 'r', x2, L21, 'g', x3, L31, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N2(x)');
xlabel('x');
ylabel('N2(x)');
grid on;

% Third subplot
subplot(2, 3, 3);
plot(x1, L12, 'r', x2, L22, 'g', x3, C3, 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N3(x)');
xlabel('x');
ylabel('N3(x)');
grid on;

% Fourth subplot
subplot(2, 3, 4);
plot(x1, C1, 'r', x2, C2, 'g', x3, L32, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N4(x)');
xlabel('x');
ylabel('N4(x)');
grid on;

% Fifth subplot
subplot(2, 3, 5);
plot(x1, L13, 'r', x2, C2, 'g', x3, C3, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N5(x)');
xlabel('x');
ylabel('N5(x)');
grid on;

% Sixth subplot
subplot(2, 3, 6);
plot(x1, C1,'r', x2, L23, 'g', x3, L33, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N6(x)');
xlabel('x');
ylabel('N6(x)');
grid on;

%% Q1-2-1
% Define the range and the number of points for plotting
x1 = linspace(0, 1, 100);
x2 = linspace(1, 2, 100);
x3 = linspace(2, 3, 100);

% Calculate the Lagrange P-3 basis functions within the ranges [0, 1], [1, 2], and [2, 3]
[L10, L11, L12, L13] = lagrange_p3_basis(x1, 0, 1);
[L20, L21, L22, L23] = lagrange_p3_basis(x2, 1, 2);
[L30, L31, L32, L33] = lagrange_p3_basis(x3, 2, 3);

C1 = zeros(100);
C2 = zeros(100);
C3 = zeros(100);

% Create the subplots
figure;

% First subplot
subplot(2, 2, 1);
plot(x1, L13, 'r', x2, L20, 'g', x3, C3, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N4(x)');
xlabel('x');
ylabel('N4(x)');
grid on;

% Second subplot
subplot(2, 2, 2);
plot(x1, C1, 'r', x2, L21, 'g', x3, C3, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N5(x)');
xlabel('x');
ylabel('N5(x)');
grid on;

% Third subplot
subplot(2, 2, 3);
plot(x1, C1, 'r', x2, L22, 'g', x3, C3, 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N6(x)');
xlabel('x');
ylabel('N6(x)');
grid on;

% Fourth subplot
subplot(2, 2, 4);
plot(x1, C1, 'r', x2, L23, 'g', x3, L30, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N7(x)');
xlabel('x');
ylabel('N7(x)');
grid on;

%% Q1-2-2

% Define the range and the number of points for plotting
x1 = linspace(0, 1, 100);
x2 = linspace(1, 2, 100);
x3 = linspace(2, 3, 100);

% Calculate the Lagrange P-3 basis functions within the ranges [0, 1], [1, 2], and [2, 3]
[L11, L12, L13, L14, L15] = lagrange_p4_basis(x1, 0, 1);
[L21, L22, L23, L24, L25] = lagrange_p4_basis(x2, 1, 2);
[L31, L32, L33, L34, L35] = lagrange_p4_basis(x3, 2, 3);

C1 = zeros(100);
C2 = zeros(100);
C3 = zeros(100);


% Create the subplots
figure;

% First subplot
subplot(2, 3, 1);
plot(x1, L15, 'r', x2, L21, 'g', x3, C3, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N5(x)');
xlabel('x');
ylabel('N5(x)');
grid on;

% Second subplot
subplot(2, 3, 2);
plot(x1, C1, 'r', x2, L22, 'g', x3, C3, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
title('N6(x)');
xlabel('x');
ylabel('N6(x)');
grid on;

% Third subplot
subplot(2, 3, 3);
plot(x1, C1, 'r', x2, L23, 'g', x3, C3, 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N7(x)');
xlabel('x');
ylabel('N7(x)');
grid on;

% Fourth subplot
subplot(2, 3, 4);
plot(x1, C1, 'r', x2, L24, 'g', x3, C3, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N8(x)');
xlabel('x');
ylabel('N8(x)');
grid on;

subplot(2, 3, 5);
plot(x1, C1, 'r', x2, L25, 'g', x3, L31, 'b', 'LineWidth', 2);
ylim([-0.4, 1.1]);
xlim([0, 3]);
title('N9(x)');
xlabel('x');
ylabel('N9(x)');
grid on;

%% q1-4

% Create the x-axis values
x = linspace(0, 3, 300);

% Define the range and the number of points for plotting
x1 = linspace(0, 1, 100);
x2 = linspace(1, 2, 100);
x3 = linspace(2, 3, 100);
x_list = [x1; x2; x3];
end_point_list = [0, 1; 1, 2; 2, 3];

% Calculate the Lagrange P-3 basis functions within the ranges [0, 1], [1, 2], and [2, 3]

[L10, L11, L12, L13] = lagrange_p3_basis(x1, 0, 1);
[L20, L21, L22, L23] = lagrange_p3_basis(x2, 1, 2);
[L30, L31, L32, L33] = lagrange_p3_basis(x3, 2, 3);

shape_func_list = {[L10, L11, L12, L13], [L20, L21, L22, L23], [L30, L31, L32, L33]};
% Create an anonymous function for f(x) = sin(n*x)
f = @(n, x) sin(n * x);
%lg_map = [ 1, 4, 7; 2, 5, 8; 3, 6, 9; 4, 7, 10];
elem_size = 3;

% Create the x-axis values
x = linspace(0, 3, 300);

% Define the values of n
n_values = [2, 4, 6];

% Define the colors for each plot
colors = ['r', 'g', 'b'];

% Loop through the values of n
for ii = 1:length(n_values)
    sol_list = zeros(size(x));
    for i = 1:elem_size
        % loop through elements
        end_point = end_point_list(i, :);
        nodes = linspace(end_point(1), end_point(2), 4);
    
        for j = 1:length(nodes)
            % loop through nodes
            f_value = f(n_values(ii), nodes(j)); % local degree of freedom
            sol_list(1, 1 + (i-1)* 100:i*100 ) = sol_list(1, 1 + (i-1)* 100:i*100 ) + f_value * shape_func_list{i}(1, (j-1)*100 + 1: j*100);
        end
    end

    
    % Create a figure
    figure(ii);
    % Calculate f(x) = sin(n*x) for the current n value
    f_n = f(n_values(ii), x);
    hold on;
    % Plot the function
    plot(x, f_n, colors(1), 'LineWidth', 2);
    plot(x, sol_list, colors(2), 'LineWidth', 2);
    plot(x, f_n - sol_list, colors(3), 'LineWidth', 2);

    % Set plot title and labels
    formattedTitle = sprintf('f(x) = sin(%dx)', n_values(ii));
    title(formattedTitle);
    xlabel('x');
    ylabel('f(x)');
    ylim([ -1, 1]);
    % Add a legend
    legend('f(x)', 'If(x)', 'f(x) - If(x)');
    % Release the hold on the current plot
    hold off;
end










