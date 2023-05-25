% hexahedral shape function visualisation in unit cube
clear all
close all
clc

% node coordinates
node_coordinates = [-1 -1 -1;
                     1 -1 -1;
                     1  1 -1;
                    -1  1 -1;
                    -1 -1  1;
                     1 -1  1;
                     1  1  1;
                    -1  1  1];

% generate xi, eta, zeta values
[xi, eta, zeta] = meshgrid(linspace(-1,1,20));

% calculate and plot shape functions
for i = 1:8
    xi_i = node_coordinates(i,1);
    eta_i = node_coordinates(i,2);
    zeta_i = node_coordinates(i,3);
    
    % shape function
    N_i = 1/8 * (1 + xi.*xi_i) .* (1 + eta.*eta_i) .* (1 + zeta.*zeta_i);
    
    % plot
    X = xi;
    Y = eta;
    Z = zeta;
    figure(i)
    % Plot slices on the boundaries of the cube
    hold on;
    surf(squeeze(X(:,:,1)), squeeze(Y(:,:,1)), squeeze(Z(:,:,1)), squeeze(N_i(:,:,1)));
    surf(squeeze(X(:,:,end)), squeeze(Y(:,:,end)), squeeze(Z(:,:,end)), squeeze(N_i(:,:,end)));
    surf(squeeze(X(:,1,:)), squeeze(Y(:,1,:)), squeeze(Z(:,1,:)), squeeze(N_i(:,1,:)));
    surf(squeeze(X(:,end,:)), squeeze(Y(:,end,:)), squeeze(Z(:,end,:)), squeeze(N_i(:,end,:)));
    surf(squeeze(X(1,:,:)), squeeze(Y(1,:,:)), squeeze(Z(1,:,:)), squeeze(N_i(1,:,:)));
    surf(squeeze(X(end,:,:)), squeeze(Y(end,:,:)), squeeze(Z(end,:,:)), squeeze(N_i(end,:,:)));

    % Configure plot
    colormap jet;
    colorbar;
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('N_i on the surface of a unit cube');

    hold off;

end

