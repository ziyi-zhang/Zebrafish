% basis plot

%knots = [1 1 1 1 2 3 4 5 6 6 6 6];
% knots = [1 1 1 1 1+5/3 1+10/3 6 6 6 6];
% knots = [1 1 1 1 1+9/7 1+18/7 1+27/7 1+36/7 1+45/7 1+54/7  10 10 10 10]; % cubic
knots = [1 1 1 1+9/7 1+18/7 1+27/7 1+36/7 1+45/7 1+54/7  10 10 10];  % quadratic

figure
hold on
grid on
degree = 3; % 3 for quadratic & 4 for cubic
for i = 1:length(knots)-degree
    b = spmak(knots(i:i+degree), 1);
    xx=1:0.01:max(knots); yy=fnval(b, xx);
    plot(xx, yy);
end
