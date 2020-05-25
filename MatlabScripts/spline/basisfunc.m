function [res] = basisfunc(t)

    if (t>1)
        res = -1/6 * t^3 + t^2 - 2*t + 4/3;
        return;
    end
    if (t>0)
        res = 0.5*t^3 - t^2 + 2/3;
        return;
    end
    if (t>-1)
        res = -0.5*t^3 - t^2 + 2/3;
        return;
    end
    % if t>-2
    res = 1/6 * t^3 + t^2 + 2*t + 4/3;
end
