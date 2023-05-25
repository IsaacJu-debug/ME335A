clear all
close all
clc

% define the vertices of a regular pentagon inscribed in a unit circle
theta = [0:4]'*2*pi/5;
vertices = [cos(theta), sin(theta)];

% generate a grid in Cartesian coordinates
[x, y] = meshgrid(linspace(-1, 1, 100));

% transform Cartesian coordinates to polar coordinates
[theta, rho] = cart2pol(x, y);

% map polar coordinates to [0, 1]
theta = mod(theta, 2*pi)/(2*pi);
rho = min(1, rho);

% calculate the shape functions
N = zeros(size(x, 1), size(x, 2), 5);
for i = 1:5
    N(:, :, i) = rho.*(0.5 + abs(mod(theta - (i-1)/5, 1) - 0.5));
end

% plot the shape functions
for i = 1:5
    figure
    surf(x, y, N(:, :, i));
    shading interp;
    colorbar;
    title(['Shape Function N' num2str(i)]);
    xlabel('x');
    ylabel('y');
    hold on;
    % draw the pentagon
    plot(vertices(:, 1), vertices(:, 2), 'k');
    plot([vertices(end, 1), vertices(1, 1)], [vertices(end, 2), vertices(1, 2)], 'k');
    hold off;
    axis equal;
end
