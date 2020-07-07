% basis function calculator

t = [1 1 1 1 2 3 4 5 6 7];  % t1, t2, t3 ...

basis = cell(5, 4);  % i, deg
syms x b11 b21 b31 b41 b51
basis{1, 1} = b11;
basis{2, 1} = b21;
basis{3, 1} = b31;
basis{4, 1} = b41;
basis{5, 1} = b51;

for deg = 2:4
    for i = 1:5-deg
        
        basis{i, deg} = (x-t(i)) / (t(i+deg-1)-t(i)) * basis{i, deg-1} + ...
                        (t(i+deg)-x) / (t(i+deg)-t(i+1)) * basis{i+1, deg-1};
    end
end

% temp = coeffs(basis{1, 4}, b41); expand(subs(temp(2), x, y+3))


