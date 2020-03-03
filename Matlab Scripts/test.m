% apply optimization on global real image

xCollect = [];
fvalCollect = [];
exitflagCollect = [];
for i = 300:20:700
    for j = 300:20:700
    
        p_ = p(i:i+29, j:j+29, :);
        f_ = @(arr)fun(double(p_), arr);
        options = optimoptions(@particleswarm, 'Display', 'off');
        [x, fval, exitflag, output] = particleswarm(f_, 5, [1 1 4 1 1], [29 29 36 6 15], options);
        xCollect = [xCollect; x]; %#ok
        fvalCollect = [fvalCollect; fval]; %#ok
        exitflagCollect = [exitflagCollect; exitflag]; %#ok
    
        fprintf('%d, %d: f:%f\n', (i-300)/20, (j-300)/20, fval);
    end
end


imshowTiff(p(300:700, 300:700, :));
hold on
count = 0;
for i = 300:20:700
    for j = 300:20:700
    
        count = count + 1;
        scatter(xCollect(count, 2)+j-300, xCollect(count, 1)+i-300, 5, 'filled');
        %viscircles([xCollect(count, 2)+j-300, xCollect(count, 1)+i-300], xCollect(count, 4), 'LineWidth', 1);
    end
end
