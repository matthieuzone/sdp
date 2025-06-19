function res = et(x, y)
    if islogical(x) && islogical(y)
        res = x && y;
    elseif islogical(x) && ~islogical(y)
        if x == true
            res = y;
        else
            res = false;
        end
    elseif ~islogical(x) && islogical(y)
        if y == true
            res = x;
        else
            res = false;
        end
    else
        res = x & y;
    end