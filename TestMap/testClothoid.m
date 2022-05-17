%testClothoid
%25/01/2020

% Parameters
a = 2;
n = 1000;
tInit = 0;
tEnd = 2;
dt = tEnd/n;

% Method 1
l = linspace(tInit,tEnd,n);
f1 = @(s)1/a*cos(s.^2);
f2 = @(s)1/a*sin(s.^2);
x = @(t)integral(f1,0,t);
y = @(t)integral(f2,0,t);
p = arrayfun(x,l);
q = arrayfun(y,l);

% Method 2
X(1) = 0;
Y(1) = 0;
t = tInit;
for i=2:n
    dx = 1/a*cos(l(i)^2) * dt;
    dy = 1/a*sin(l(i)^2) * dt;
    X(i) = (X(i-1) + dx);
    Y(i) = (Y(i-1) + dy);
    t = t + dt;
end

% Plots
figure(1);
plot(p,q,'b.-');
hold on;
plot(X,Y,'r.-');
legend('Method 1', 'Method 2');
title('Clothoid');

figure(2);
plot(l, (abs(p)-abs(X)),'b.-');
title('error X');

figure(3);
plot(l, (abs(q)-abs(Y)),'b.-');
title('error Y');
