function [res] = fun(mat, cylArray, z)

    cyl.x = cylArray(1);
    cyl.y = cylArray(2);
    cyl.z = z;  % no longer optimize 'z'
    cyl.r = cylArray(3);
    cyl.h = 4;  % fix height as 4
    
      % Do this once before calling this function:
      % mat = double(mat);
      % mat = mat ./ max(mat, [], 'all');  % normalize to 0-1
    [sampleCyl, samplePeri, weight] = SampleCylinder(cyl, mat, false, 'equadistant');
    res = EvaluateCylinder(mat, sampleCyl, samplePeri, weight);
end
