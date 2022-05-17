%testMap
%26/04/2019

close all;
clear;
clc;
 
% Map points
% Segmet 1
p1 = [0; 0];
p2 = [0; 25];
p3 = [0; 75];
p4 = [0; 100];

% Segmet 2
p5 = p4;
p6 = [25;  100];
p7 = [50;  100];
p8 = [100; 100];

segment1 = Bezier3('nPoints',500,'p1',p1,'p2',p2,'p3',p3,'p4',p4);
segment2 = Bezier3('nPoints',500,'p1',p5,'p2',p6,'p3',p7,'p4',p8);
segment1 = segment1.getPoints(segment1);
segment2 = segment2.getPoints(segment2);
segment1.plotBezier(segment1);
hold on;
segment2.plotBezier(segment2);


Map.x = [segment1.points(1,1:end-1),segment2.points(1,:)];
Map.y = [segment1.points(2,1:end-1),segment2.points(2,:)];


clear distPoints  p1 p2 p3 p4 p5 p6 p7 p8 segment1 segment2
save Map.mat

clear Map;
 
