function [L0, L1, L2] = lagrange_p2_basis(x, x0, x1)
    % Define the reference nodes for Lagrange P-2 element
    nodes = linspace(x0, x1, 3);

    % Calculate the Lagrange basis functions
    L0 = (x - nodes(2)) .* (x - nodes(3)) / ...
        ((nodes(1) - nodes(2)) * (nodes(1) - nodes(3)));

    L1 = (x - nodes(1)) .* (x - nodes(3)) / ...
        ((nodes(2) - nodes(1)) * (nodes(2) - nodes(3)));

    L2 = (x - nodes(1)) .* (x - nodes(2)) / ...
        ((nodes(3) - nodes(1)) * (nodes(3) - nodes(2)));
end
