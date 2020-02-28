function [res] = fun(mat, cylArray)

    cyl.x = cylArray(1);
    cyl.y = cylArray(2);
    cyl.z = cylArray(3);
    cyl.r = cylArray(4);
    cyl.h = cylArray(5);
    [sampleCyl, samplePeri] = SampleCylinder(cyl);
    res = EvaluateCylinder(mat, sampleCyl, samplePeri);
end

% testMat = ones(100, 100, 30);
% testMat(20:22, 30:32, 6:10) = 0;