%testPath
%24/06/2019
close all
clear
clc


path = Path('Bezier3');
path.plotPath(path,100);
t = linspace(0,1,1000);
x0 = 0;
y0 = 0;

for i = 1:length(t)
    [x(i) y(i)] = path.getPath(path,t(i),2);
    f(i)   =((x(i)-x0)^2) + ((y(i)-y0)^2) - 10^2;
end

plot(x,f,'b.-');