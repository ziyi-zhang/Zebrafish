function [res] = fun(mat, cylArray)

    cyl.x = cylArray(1);
    cyl.y = cylArray(2);
    cyl.z = cylArray(3);
    cyl.r = cylArray(4);
    % cyl.h = cylArray(5);
    mat = mat ./ max(mat, [], 'all');  % normalize to 0-1
    [sampleCyl, samplePeri, weight] = SampleCylinder(cyl, mat, false, 'equadistant');
    res = EvaluateCylinder(mat, sampleCyl, samplePeri, weight);
end
