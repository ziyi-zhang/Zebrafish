
temp = temp';  % now 57*3
fprintf("\n");
count = 0;
for i = 1:size(temp, 1)
    for j = 1:size(temp, 2)
    
        fprintf("%f, ", temp(i, j));
        count = count + 1;
        if (mod(count, 500) == 0)
            fprintf("\n");
        end
    end
    %fprintf("\n");
end
fprintf("\n");