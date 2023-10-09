% quadrilateral shape function visualisation in unit square
clear all
close all
clc

% node coordinates
node_coordinates = [-1 -1;
                     1 -1;
                     1  1;
                    -1  1];

% generate xi and eta values
[xi, eta] = meshgrid(linspace(-1,1,20));

% calculate and plot shape functions
for i = 1:4
    xi_i = node_coordinates(i,1);
    eta_i = node_coordinates(i,2);
    
    % shape function
    N_i = 1/4 * (1 + xi*xi_i) .* (1 + eta*eta_i);
    
    % plot
    figure
    surf(xi, eta, N_i)
    title(['Shape Function N' num2str(i)])
    xlabel('xi')
    ylabel('eta')
    zlabel(['N' num2str(i)])
    axis equal
    axis([-1 1 -1 1 0 1])
    colorbar
    grid on
end
