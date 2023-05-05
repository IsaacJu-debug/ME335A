function [L0, L1, L2, L3] = lagrange_p3_basis(x, x0, x1)
    % Define the reference nodes for Lagrange P-3 element
    nodes = linspace(x0, x1, 4);

    % Calculate the Lagrange basis functions
    L0 = (x - nodes(2)) .* (x - nodes(3)) .* (x - nodes(4)) / ...
        ((nodes(1) - nodes(2)) * (nodes(1) - nodes(3)) * (nodes(1) - nodes(4)));

    L1 = (x - nodes(1)) .* (x - nodes(3)) .* (x - nodes(4)) / ...
        ((nodes(2) - nodes(1)) * (nodes(2) - nodes(3)) * (nodes(2) - nodes(4)));

    L2 = (x - nodes(1)) .* (x - nodes(2)) .* (x - nodes(4)) / ...
        ((nodes(3) - nodes(1)) * (nodes(3) - nodes(2)) * (nodes(3) - nodes(4)));

    L3 = (x - nodes(1)) .* (x - nodes(2)) .* (x - nodes(3)) / ...
        ((nodes(4) - nodes(1)) * (nodes(4) - nodes(2)) * (nodes(4) - nodes(3)));
end
