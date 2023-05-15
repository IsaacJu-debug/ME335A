function [N1, N2, N3] = p1_triangle(x, y, x1, y1, x2, y2, x3, y3)
    % Calculate the shape functions for a triangular element
    % Input:
    % x, y : Cartesian coordinates of the point where the shape function is to be evaluated
    % x1, y1, x2, y2, x3, y3 : Cartesian coordinates of the triangle's vertices
    %
    % Output:
    % N1, N2, N3 : Values of shape functions at (x, y)

    % Calculate the shape functions
    N1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
    N2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
    N3 = 1 - N1 - N2;
end