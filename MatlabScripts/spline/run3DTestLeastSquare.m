
f = @(x, y, z) x.^2 + x.*z - y.*x;

[x, y, z] = meshgrid(1:4, 1:4, 1:4);
feval = f(x, y, z);

% Least square fitting
N_x = size(z, 1);
N_y = size(z, 2);
N_z = size(z, 3);
N = N_x * N_y * N_z;
num_x = 4;
num_y = 4;
num_z = 4;
num = num_x * num_y * num_z;
gap_x = (N_x-1) / (num_x-1);
gap_y = (N_y-1) / (num_y-1);
gap_z = (N_z-1) / (num_z-1);
centers_x = (0:num_x-1).*gap_x + 1; 
centers_y = (0:num_y-1).*gap_y + 1;
centers_z = (0:num_z-1).*gap_z + 1;
A = zeros(N, num);
for jz = 1:num_z
    for jx = 1:num_x
        for jy = 1:num_y
            j = (jz-1)*num_y*num_x + (jx-1)*num_y + jy;
            center_x = centers_x(jx);
            center_y = centers_y(jy);
            center_z = centers_z(jz);
            for iz = 1:N_z
                for iy = 1:N_y
                    for ix = 1:N_x
                        i = (iz-1)*N_y*N_x + (iy-1)*N_x + ix;

                        A(i, j) = basisfunc((ix-center_x)/gap_x) *...
                                  basisfunc((iy-center_y)/gap_y) *...
                                  basisfunc((iz-center_z)/gap_z);
                              
                              
                        % a special fix to the border issue
                        %{
                        count = 0;
                        if (ix==1 || ix==num_x), count = count+1;end
                        if (iy==1 || iy==num_y), count = count+1;end
                        if (iz==1 || iz==num_z), count = count+1;end
                        if (count == 1)
                            A(i, j) = A(i, j) / 0.8333;
                        elseif (count == 2)
                            A(i, j) = A(i, j) / 0.6944;
                        elseif (count == 3)
                            A(i, j) = A(i, j) / 0.5787;
                        end
                        %}
                    end
                end
            end
        end
    end
end
inputPts = feval(:);
controlPts = linsolve(A'*A, A'*inputPts);

% test
for i = 1:10

    % query = randi(10, 1, 3) + 5;  % intersections
    query = rand(1, 3) .* 10 + 5;
    res1 = interp3DLeastSquare(centers_x, centers_y, centers_z, gap_x, gap_y, gap_z, controlPts, query);
    res2 = interp3(feval, query(2), query(1), query(3), 'spline');  % pass in the inputPts, not controlPts
    roundquery = round(query);
    avg = feval(roundquery(1)-1:roundquery(1)+1, roundquery(2)-1:roundquery(2)+1, roundquery(3)-1:roundquery(3)+1);
    avg = mean(avg(:));
    fprintf("mean=%10.3f | b-spline(fitting)=%10.3f | Interp3=%10.3f\n", avg, res1, res2);
end
