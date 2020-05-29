function [res] = interp3DLeastSquare(centers_x, centers_y, centers_z, gap_x, gap_y, gap_z, controlPts, query)

    num_x = length(centers_x);
    num_y = length(centers_y);
    num_z = length(centers_z);
    idx_x = floor( (query(1)-1)/gap_x );  % lowest idx of the 4*4*4 cube
    idx_y = floor( (query(2)-1)/gap_y );
    idx_z = floor( (query(3)-1)/gap_z );
    
    res = 0;
    for iz = 0:3
        for ix = 0:3
            for iy = 0:3
            
                x = idx_x + ix;
                y = idx_y + iy;
                z = idx_z + iz;
                arrayIdx = num_x*num_y*(z-1) + num_y*(x-1)+y;
                
                xcoef = basisfunc( (query(1)-centers_x(x)) / gap_x );
                ycoef = basisfunc( (query(2)-centers_y(y)) / gap_y );
                zcoef = basisfunc( (query(3)-centers_z(z)) / gap_z );
                
                res = res + controlPts(arrayIdx) * xcoef * ycoef * zcoef;
            end
        end
    end
end

% vector query?
% matrix BxyzIdx?
