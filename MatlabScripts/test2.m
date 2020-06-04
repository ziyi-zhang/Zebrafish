
temp = temp';  % now 57*3
fprintf("\n");
for i = 1:size(temp, 1)
    for j = 1:size(temp, 2)
    
        fprintf("%f, ", temp(i, j));
    end
end
fprintf("\n");