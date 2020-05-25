function [res] = interp3D(image, query)

    truncQ = floor(query);
    baseIdx = truncQ - 1;
    
    z1layer = image(baseIdx(1):baseIdx(1)+3, baseIdx(2):baseIdx(2)+3, baseIdx(3)  );
    z2layer = image(baseIdx(1):baseIdx(1)+3, baseIdx(2):baseIdx(2)+3, baseIdx(3)+1);
    z3layer = image(baseIdx(1):baseIdx(1)+3, baseIdx(2):baseIdx(2)+3, baseIdx(3)+2);
    z4layer = image(baseIdx(1):baseIdx(1)+3, baseIdx(2):baseIdx(2)+3, baseIdx(3)+3);
    yArray = [z1layer(:)', z2layer(:)', z3layer(:)', z4layer(:)'];
    
    % Bx
    BxIdx = repmat([0 1 2 3], 1, 16);
    t1 = query(1) - truncQ(1) + 1;
    t2 = query(1) - truncQ(1);
    t3 = query(1) - truncQ(1) - 1;
    t4 = query(1) - truncQ(1) - 2;
    BxArray = [-1/6*t1^3+t1^2-2*t1+4/3, ...
               0.5*t2^3-t2^2+2/3, ...
               -0.5*t3^3-t3^2+2/3, ...
               1/6*t4^3+t4^2+2*t4+4/3];
    BxArray = BxArray(BxIdx + 1);
    
    % By
    ByIdx = repmat([0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3], 1, 4);
    t1 = query(2) - truncQ(2) + 1;
    t2 = query(2) - truncQ(2);
    t3 = query(2) - truncQ(2) - 1;
    t4 = query(2) - truncQ(2) - 2;
    ByArray = [-1/6*t1^3+t1^2-2*t1+4/3, ...
               0.5*t2^3-t2^2+2/3, ...
               -0.5*t3^3-t3^2+2/3, ...
               1/6*t4^3+t4^2+2*t4+4/3];
    ByArray = ByArray(ByIdx + 1);
    
    % Bz
    temp = ones(1, 16);
    BzIdx = [0.*temp, temp, 2.*temp, 3.*temp];
    t1 = query(3) - truncQ(3) + 1;
    t2 = query(3) - truncQ(3);
    t3 = query(3) - truncQ(3) - 1;
    t4 = query(3) - truncQ(3) - 2;
    BzArray = [-1/6*t1^3+t1^2-2*t1+4/3, ...
               0.5*t2^3-t2^2+2/3, ...
               -0.5*t3^3-t3^2+2/3, ...
               1/6*t4^3+t4^2+2*t4+4/3];
    BzArray = BzArray(BzIdx + 1);
    
    res = sum(yArray .* BxArray .* ByArray .* BzArray);
end

% vector query?
% matrix BxyzIdx?
