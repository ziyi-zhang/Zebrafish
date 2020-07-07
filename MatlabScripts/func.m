function [res] = func(x)

    res = x - 1;
    
    mask = x < 0;
    res(mask) = res(mask) * (-1);
end
