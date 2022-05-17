%testMap
%26/04/2019

close all;
clear;
clc;
 
% Map points
% Segmet 1
p1 = [0; 0];
p2 = [0; 500];
p3 = [500; 0];
p4 = [500; 500];

segment1 = Bezier3('nPoints',500,'p1',p1,'p2',p2,'p3',p3,'p4',p4);

segment1 = segment1.getPoints(segment1);

segment1.plotBezier(segment1);

Map.x = [segment1.points(1,1:end)];
Map.y = [segment1.points(2,1:end)];


clear distPoints  p1 p2 p3 p4 segment1
save Map.mat

clear Map;
 
