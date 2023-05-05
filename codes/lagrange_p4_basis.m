function [L0, L1, L2, L3, L4] = lagrange_p4_basis(x, x0, x1)
    % Define the reference nodes for Lagrange P-4 element
    nodes = linspace(x0, x1, 5);

    % Calculate the Lagrange basis functions
    L0 = (x - nodes(2)) .* (x - nodes(3)) .* (x - nodes(4)) .* (x - nodes(5)) / ...
        ((nodes(1) - nodes(2)) * (nodes(1) - nodes(3)) * (nodes(1) - nodes(4)) * (nodes(1) - nodes(5)));

    L1 = (x - nodes(1)) .* (x - nodes(3)) .* (x - nodes(4)) .* (x - nodes(5)) / ...
        ((nodes(2) - nodes(1)) * (nodes(2) - nodes(3)) * (nodes(2) - nodes(4)) * (nodes(2) - nodes(5)));

    L2 = (x - nodes(1)) .* (x - nodes(2)) .* (x - nodes(4)) .* (x - nodes(5)) / ...
        ((nodes(3) - nodes(1)) * (nodes(3) - nodes(2)) * (nodes(3) - nodes(4)) * (nodes(3) - nodes(5)));

    L3 = (x - nodes(1)) .* (x - nodes(2)) .* (x - nodes(3)) .* (x - nodes(5)) / ...
        ((nodes(4) - nodes(1)) * (nodes(4) - nodes(2)) * (nodes(4) - nodes(3)) * (nodes(4) - nodes(5)));

    L4 = (x - nodes(1)) .* (x - nodes(2)) .* (x - nodes(3)) .* (x - nodes(4)) / ...
        ((nodes(5) - nodes(1)) * (nodes(5) - nodes(2)) * (nodes(5) - nodes(3)) * (nodes(5) - nodes(4)));
end
