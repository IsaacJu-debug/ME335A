%% q1-p3
clear
syms x_1 x_2

f1 = sin(pi * (x_1 ^ 2 + x_2 ^ 2));
f2 = cos(pi/2 * (x_1 ^ 2 + x_2 ^ 2));

% Compute the derivative of N1 with respect to x_1
df1_dx1 = diff(f1, x_1);
disp(df1_dx1); 
% Compute the derivative of N1 with respect to x_2
df1_dx2 = diff(f1, x_2);
disp(df1_dx2);

% Compute the derivative of N2 with respect to x_1
df2_dx1 = diff(f2, x_1);
disp(df2_dx1); 
% Compute the derivative of N2 with respect to x_2
df2_dx2 = diff(f2, x_2);
disp(df2_dx2);

% compute the bilinear term
df1_dx = [df1_dx1, df1_dx2];
disp(df1_dx);

df2_dx = [df2_dx1, df2_dx2];
disp(df2_dx);

% Compute the dot product of df_dx and df_dx
df1_dot = dot(df1_dx, df1_dx);
df2_dot = dot(df2_dx, df2_dx);
df1_df2_dot = dot(df1_dx, df2_dx);

% Convert to polar coordinates for integration over the disk
syms r theta
df1_dot_polar = subs(df1_dot, {x_1, x_2}, {r*cos(theta), r*sin(theta)});
df1_dot_polar = simplify(df1_dot_polar);

% Compute the integral over the disk
I_N1_N1 = int(int(r*df1_dot_polar, r, 0, 1), theta, 0, 2*pi);
%I_N1_N1 = int(int(df1_dot_polar, r, 0, 1), theta, 0, 2*pi);

% Multiply the result by 1/2
I_N1_N1 = I_N1_N1 / 2;
disp(I_N1_N1);
disp(double(I_N1_N1));

%% N1_N2
df1_df2_dot_polar = subs(df1_df2_dot, {x_1, x_2}, {r*cos(theta), r*sin(theta)});
df1_df2_dot_polar = simplify(df1_df2_dot_polar);

% Compute the integral over the disk
I_N1_N2 = int(int(r*df1_df2_dot_polar, r, 0, 1), theta, 0, 2*pi);

% Multiply the result by 1/2
I_N1_N2 = I_N1_N2 / 2;
disp(double(I_N1_N2));

%% N2_N2
df2_dot_polar = subs(df2_dot, {x_1, x_2}, {r*cos(theta), r*sin(theta)});
df2_dot_polar = simplify(df2_dot_polar);

% Compute the integral over the disk
I_N2_N2 = int(int(r*df2_dot_polar, r, 0, 1), theta, 0, 2*pi);

% Multiply the result by 1/2
I_N2_N2 = I_N2_N2 / 2;
disp(I_N2_N2);
disp(double(I_N2_N2));
%% l(N1) = int(2/R^2 * N1)_Omega - int(N1/R)_Gamma_h
f1_polar = subs(f1, {x_1, x_2}, {r*cos(theta), r*sin(theta)});
% Compute the integral over the disk
I_N1_1 = 2 * int(int(r*f1_polar, r, 0, 1), theta, 0, 2*pi);

% The radius r is 1 on the unit circle
N1_polar = 2.0 *subs(f1_polar, r, 1);

% Compute the integral over the arc of the unit circle from (pi/2, 2*pi)
I_N1_2 = -int(N1_polar, theta, pi/2, 2*pi);

I_N1 = I_N1_2 + I_N1_1;
disp(I_N1);

%% l(N2)
f2_polar = subs(f2, {x_1, x_2}, {r*cos(theta), r*sin(theta)});
% Compute the integral over the disk
I_N2_1 = 2 * int(int(r*f2_polar, r, 0, 1), theta, 0, 2*pi);

% The radius r is 1 on the unit circle
N2_polar =  2.0 *subs(f2_polar, r, 1);

% Compute the integral over the arc of the unit circle from (pi/2, 2*pi)
I_N2_2 = -int(N2_polar, theta, pi/2, 2*pi);

I_N2 = I_N2_2 + I_N2_1;
disp(I_N2);

%% q2 solve the system
K = [ pi^3/2, 20 * pi/ 9; 20 * pi/ 9, (pi^3 + 4*pi)/8];
F = [4, 4]';
U = K\F;
disp(U);
%U = [-0.1720, 0.9548]';
% Define the symbolic variables and functions
syms x_1 x_2
N1 = sin(pi * (x_1 ^ 2 + x_2 ^ 2));
N2 = cos(pi/2 * (x_1 ^ 2 + x_2 ^ 2));

% Define u_h
u_h = U(1)*N1 + U(2)*N2;

% Convert to a function handle for numerical computation
u_h_func = matlabFunction(u_h);

% Define a grid of x1 and x2 values
[x1, x2] = meshgrid(linspace(-1, 1, 1000), linspace(-1, 1, 1000));

% Compute the function values at each grid point
u_h_values = zeros(size(x1));
u_values = zeros(size(x1));
for i = 1:numel(x1)
    if x1(i)^2 + x2(i)^2 <= 1  % only compute values inside the unit circle
        u_h_values(i) = u_h_func(x1(i), x2(i));
        u_values(i) = 1 - x1(i)^2 - x2(i)^2;
    else
        u_h_values(i) = 0.0;  % assign NaN outside the unit circle
        u_values(i) = 0.0;
    end
end

% Create a surface plot of the function over the unit circle
figure(1);
surf(x1, x2, u_h_values, 'EdgeColor', 'none');
xlabel('x1');
ylabel('x2');
zlabel('u_h');
title('FE solution over the unit circle');
colorbar;

figure(2);
surf(x1, x2, u_values, 'EdgeColor', 'none');
xlabel('x1');
ylabel('x2');
zlabel('u_h');
title('Exact solution over the unit circle');
colorbar;

figure(3);
surf(x1, x2, u_values - u_h_values, 'EdgeColor', 'none');
xlabel('x1');
ylabel('x2');
zlabel('u_h');
title('Absolute error the unit circle');
colorbar;



