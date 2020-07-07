% % basis function calculator record

x21 = -2:0.01:-1;
x10 = -1:0.01:0;
x01 = 0:0.01:1;
x12 = 1:0.01:2;
x23 = 2:0.01:3;
x34 = 3:0.01:4;


% 1
%{
% cubic 111
temp2 = (x/2 - 3/2)*((x/2 - 3/2)*(b31*(x - 1) - b41*(x - 3)) - (0*b21*(x - 1) - b31*(x - 2))*(x - 1)) + (0*(x - 1)*(0*b11*(x - 1) - 0*b21*(x - 1)) - (0*b21*(x - 1) - b31*(x - 2))*(x - 2))*(x - 1);

temp = coeffs(temp2, b31); b31coeff = expand(subs(temp(2), x, y+1));
temp = coeffs(temp2, b41); b41coeff = expand(subs(temp(2), x, y+1));
b31y = double(subs(b31coeff, y, x01));
b41y = double(subs(b41coeff, y, x12));
% plot
figure
hold on
grid on
plot(x01, b31y)
plot(x12, b41y)
%}

% 2
%
% cubic 11
temp2 = - (x/3 - 4/3)*((x/2 - 1/2)*(b31*(x - 1) - b41*(x - 3)) - (x/2 - 2)*(b41*(x - 2) - b51*(x - 4))) - (x/2 - 1/2)*((x/2 - 3/2)*(b31*(x - 1) - b41*(x - 3)) - (0*b21*(x - 1) - b31*(x - 2))*(x - 1));

temp = coeffs(temp2, b31); b31coeff = expand(subs(temp(2), x, y+1));
temp = coeffs(temp2, b41); b41coeff = expand(subs(temp(2), x, y+1));
temp = coeffs(temp2, b51); b51coeff = expand(subs(temp(2), x, y+1));
b31y = double(subs(b31coeff, y, x01));
b41y = double(subs(b41coeff, y, x12));
b51y = double(subs(b51coeff, y, x23));
% plot
figure
hold on
grid on
plot(x01, b31y)
plot(x12, b41y)
plot(x23, b51y)
%

% 3
%{
% cubic 1234
temp2 = (x/3 - 1/3)*((x/2 - 1/2)*(b11*(x - 1) - b21*(x - 3)) - (x/2 - 2)*(b21*(x - 2) - b31*(x - 4))) - (x/3 - 5/3)*((x/2 - 1)*(b21*(x - 2) - b31*(x - 4)) - (x/2 - 5/2)*(b31*(x - 3) - b41*(x - 5)));

temp = coeffs(temp2, b11); b11coeff = expand(subs(temp(2), x, y+3));
temp = coeffs(temp2, b21); b21coeff = expand(subs(temp(2), x, y+3));
temp = coeffs(temp2, b31); b31coeff = expand(subs(temp(2), x, y+3));
temp = coeffs(temp2, b41); b41coeff = expand(subs(temp(2), x, y+3));
y21 = double(subs(b11coeff, y, x21));
y10 = double(subs(b21coeff, y, x10));
y01 = double(subs(b31coeff, y, x01));
y12 = double(subs(b41coeff, y, x12));
% plot
figure
hold on
grid on
plot(x21, y21)
plot(x10, y10)
plot(x01, y01)
plot(x12, y12)
%}


% 4
%{
% quadratic 1234
temp2 = (x/2 - 1/2)*(b11*(x - 1) - b21*(x - 3)) - (x/2 - 2)*(b21*(x - 2) - b31*(x - 4));

temp = coeffs(temp2, b11); b11coeff = expand(subs(temp(2), x, y+2));
temp = coeffs(temp2, b21); b21coeff = expand(subs(temp(2), x, y+2));
temp = coeffs(temp2, b31); b31coeff = expand(subs(temp(2), x, y+2));
y10 = double(subs(b11coeff, y, x10));
y01 = double(subs(b21coeff, y, x01));
y12 = double(subs(b31coeff, y, x12));
% plot
figure
hold on
grid on
plot(x10, y10)
plot(x01, y01)
plot(x12, y12)
%}


% 5
%{
% quadratic 11234
temp2 = (0*b11*(x - 1) - b21*(x - 2))*(x - 1) - (x/2 - 3/2)*(b21*(x - 1) - b31*(x - 3));
temp = coeffs(temp2, b21); b21coeff = expand(subs(temp(2), x, y+2));
temp = coeffs(temp2, b31); b31coeff = expand(subs(temp(2), x, y+2));
y10 = double(subs(b21coeff, y, x10));
y01 = double(subs(b31coeff, y, x01));
% plot
figure
hold on
grid on
plot(x10, y10)
plot(x01, y01)
%}

