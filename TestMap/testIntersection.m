% testIntersectionLines
% 19/11/2020

close all;
clear;
clc;

P1 = [4; 2];
P1_ =[1; 0];

P2 = [4; 5];
P2_ = [3; 0];

x1 = [P1_(1) P1(1)];
y1 = [P1_(2) P1(2)];

x2 = [P2(1) P2_(1)];
y2 = [P2(2) P2_(2)];

%fit linear polynomial
p1 = polyfit(x1,y1,1);
p2 = polyfit(x2,y2,1);
%calculate intersection
x_intersect = fzero(@(x) polyval(p1-p2,x),3);
y_intersect = polyval(p1,x_intersect);
line(x1,y1);
hold on;
line(x2,y2);
plot(x_intersect,y_intersect,'r*');

