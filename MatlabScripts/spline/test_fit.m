close all;

nsamples = 12;
nctrl = 12;
degree = 3;

range = 1;


% func = @(x,y)(x.^3+x.*y);
func = @(x, y)((x-14).^3 + (y-15).^3);
% func = @(x, y) ((x).^2 + (y).^2);

[xs, ys] = meshgrid(linspace(0,range,nsamples),linspace(0,range,nsamples));

x = xs(:);
y = ys(:);
v = func(x,y);

knots = linspace(0, range, nctrl+degree+1-degree*2);
knots = [ones(1, degree)*knots(1) knots  ones(1, degree)*knots(end)];
bases = cell(nctrl*nctrl, 1);

index = 1;
for i=1:nctrl
    knotsu = knots(i:i+degree+1);
    for j=1:nctrl
        knotsv = knots(j:j+degree+1);
        
        % fprintf("%d %d: ", i, j);
        % knotsv
        bases{index} = spmak({knotsu, knotsv},1);
        assert(bases{index}.order(1) == degree+1);
        assert(bases{index}.order(2) == degree+1);
        index = index+1;
    end
end

A = zeros(nsamples*nsamples,nctrl*nctrl)*nan;

for i=1:nctrl*nctrl
    A(:, i) = fnval(bases{i},[x,y]')';
end
A=sparse(A);

mat = A'*A;
ctrl = mat\ (A'*v);

interp = A*ctrl;

err = interp-v;

display(median(abs(err)));
    

figure;
surf(reshape(x,nsamples,nsamples),reshape(y,nsamples,nsamples),reshape(err,nsamples,nsamples));
hold on
% surf(reshape(x,nsamples,nsamples),reshape(y,nsamples,nsamples),reshape(v,nsamples,nsamples));
% surf(reshape(x,nsamples,nsamples),reshape(y,nsamples,nsamples),reshape(interp,nsamples,nsamples));

