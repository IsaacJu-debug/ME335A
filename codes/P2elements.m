function [y]=P2elements(x, x1, x2, x3)
    t = size(x);
    y = zeros([size(x)]); 
    for ii=1:t(2)
        if x(ii) >= min([x1, x2, x3]) && x(ii) <= max([x1, x2, x3])
            y(ii) = (x(ii) - x2)*(x(ii) - x3)/((x1 - x2)*(x1 - x3));
        else
            y(ii) = 0;
        end
    end
end