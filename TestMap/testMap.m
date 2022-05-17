%testMap
%19/11/2020

close all;
clear;
clc;

addpath('functions');
addpath('data');

% Circles
C1.x = 150;
C1.y = 300;
C1.r = 60;

C2.x = 700;
C2.y = 200;
C2.r = 120;

circle1.yaw = linspace(0, 2*pi, 100);
circle1.x = C1.x + (C1.r * cos(circle1.yaw));
circle1.y = C1.y + (C1.r * sin(circle1.yaw));
circle2.yaw = linspace(0, 2*pi, 200);
circle2.x = C2.x + (C2.r * cos(circle2.yaw));
circle2.y = C2.y + (C2.r * sin(circle2.yaw));

% Circle junction points
F1.x = circle1.x(40);
F1.y = circle1.y(40);
F1.yaw = circle1.yaw(40) + pi/2;
F1.u = cos(F1.yaw);
F1.v = sin(F1.yaw);

F2.x = circle1.x(60);
F2.y = circle1.y(60);
F2.yaw = circle1.yaw(60) - pi/2;
F2.u = cos(F2.yaw);
F2.v = sin(F2.yaw);

F3.x = circle2.x(40);
F3.y = circle2.y(40);
F3.yaw = circle2.yaw(40) - pi/2;
F3.u = cos(F3.yaw);
F3.v = sin(F3.yaw);

F4.x = circle2.x(150);
F4.y = circle2.y(150);
F4.yaw = circle2.yaw(150) + pi/2;
F4.u = cos(F4.yaw);
F4.v = sin(F4.yaw);

% Clothoidal segments
phi1_0 = deg2rad(190);
phi2_0 = deg2rad(150);
phi3_0 = deg2rad(20);
phi4_0 = deg2rad(320);

[segment1.x, segment1.y] = getClothoidPoints([F1.x; F1.y], [C1.x; C1.y],  1, phi1_0, F1.yaw);
[segment2.x, segment2.y] = getClothoidPoints([F2.x; F2.y], [C1.x; C1.y], -1, phi2_0, F2.yaw);
[segment3.x, segment3.y] = getClothoidPoints([F3.x; F3.y], [C2.x; C2.y], -1, phi3_0, F3.yaw);
[segment4.x, segment4.y] = getClothoidPoints([F4.x; F4.y], [C2.x; C2.y],  1, phi4_0, F4.yaw);

% Lines
Q1 = [segment1.x(1), segment1.y(1)];
Q2 = [segment2.x(1), segment2.y(1)];
Q3 = [segment3.x(1), segment3.y(1)];
Q4 = [segment4.x(1), segment4.y(1)];

[Q1_(1), Q1_(2)] = getEndLinePoint(Q1(1), Q1(2), 700, phi1_0);
[Q2_(1), Q2_(2)] = getEndLinePoint(Q2(1), Q2(2), 500, phi2_0);
[Q3_(1), Q3_(2)] = getEndLinePoint(Q3(1), Q3(2), 100, phi3_0);
[Q4_(1), Q4_(2)] = getEndLinePoint(Q4(1), Q4(2), 100, phi4_0);

line1.x1 = Q1(1);
line1.y1 = Q1(2);
line1.x2 = Q1_(1);
line1.y2 = Q1_(2);
line2.x1 = Q2(1);
line2.y1 = Q2(2);
line2.x2 = Q2_(1);
line2.y2 = Q2_(2);
line3.x1 = Q3(1);
line3.y1 = Q3(2);
line3.x2 = Q3_(1);
line3.y2 = Q3_(2);
line4.x1 = Q4(1);
line4.y1 = Q4(2);
line4.x2 = Q4_(1);
line4.y2 = Q4_(2);

% Clothoids to connect lines
V1 = getIntersectionPoint(line1, line4);
V2 = getIntersectionPoint(line2, line3);

phi1_1 = deg2rad(getAngle(Q1, V1));
phi2_1 = deg2rad(getAngle(Q4, V1));
phi3_1 = deg2rad(getAngle(V2, Q2)) + pi;
phi4_1 = deg2rad(getAngle(V2, Q3)) + pi;
[segment5, segment6] = getClothoidCircleClothoidPoints(Q1, Q4, V1, -1, 1, 50, phi1_1, phi2_1);
[segment7, segment8] = getClothoidCircleClothoidPoints(Q2, Q3, V2, 1, -1, 70, phi3_1, phi4_1);

% Straight lines
x = linspace(Q1(1),segment5.x(1), 100);
[~, y] = getEndLinePoint(Q1(1),Q1(2), x, phi1_0);
segment9.x = x(2:end-1);
segment9.y = y(2:end-1);
x = linspace(Q2(1),segment7.x(1), 100);
[~, y] = getEndLinePoint(Q2(1),Q2(2), x, phi2_0);
segment10.x = x(2:end-1);
segment10.y = y(2:end-1);
x = linspace(Q3(1),segment8.x(1), 100);
[~, y] = getEndLinePoint(Q3(1),Q3(2), x, phi3_0);
segment11.x = x(2:end-1);
segment11.y = y(2:end-1);
x = linspace(Q4(1),segment6.x(1), 100);
[~, y] = getEndLinePoint(Q4(1),Q4(2), x, phi4_0);
segment12.x = x(2:end-1);
segment12.y = y(2:end-1);

% Circular segments
segment13.x = circle1.x(40+1:60-1);
segment13.y = circle1.y(40+1:60-1);

segment14.x = [circle2.x(150+1:199) circle2.x(1:40-1)];
segment14.y = [circle2.y(150+1:199) circle2.y(1:40-1)];


% Pack map
s1.x = flip(segment1.x);
s1.y = flip(segment1.y);
s2.x = (segment9.x)';
s2.y = (segment9.y)';
s3.x = (segment5.x);
s3.y = (segment5.y);
s4.x = flip(segment6.x);
s4.y = flip(segment6.y);
s5.x = flip(segment12.x)';
s5.y = flip(segment12.y)';
s6.x = (segment4.x);
s6.y = (segment4.y);
s7.x = (segment14.x)';
s7.y = (segment14.y)';
s8.x = flip(segment3.x);
s8.y = flip(segment3.y);
s9.x = (segment11.x)';
s9.y = (segment11.y)';
s10.x = (segment8.x);
s10.y = (segment8.y);
s11.x = flip(segment7.x);
s11.y = flip(segment7.y);
s12.x = flip(segment10.x)';
s12.y = flip(segment10.y)';
s13.x = (segment2.x);
s13.y = (segment2.y);
s14.x = flip(segment13.x)';
s14.y = flip(segment13.y)';

map.x = vertcat(s1.x, s2.x, s3.x, s4.x, s5.x, s6.x, s7.x, s8.x, s9.x, s10.x, s11.x, s12.x, s13.x, s14.x);
map.y = vertcat(s1.y, s2.y, s3.y, s4.y, s5.y, s6.y, s7.y, s8.y, s9.y, s10.y, s11.y, s12.y, s13.y, s14.y);

Map.x = map.x';
Map.y = map.y';

save('data\loopMap.mat', 'Map');


figure(1);
hold on;
grid on;
% xlim([0, 800]);
% ylim([0, 600]);
xlabel('x [m]');
ylabel('y [m]');
% quiver(P3.x, P3.y, P3.u, P3.v, 'AutoScaleFactor', 20);
% quiver(P4.x, P4.y, P4.u, P4.v, 'AutoScaleFactor', 20);
plot(circle1.x, circle1.y, 'k.-');
plot(circle2.x, circle2.y, 'k.-');
% quiver(F1.x, F1.y, F1.u, F1.v, 'AutoScaleFactor', 20);
% quiver(F2.x, F2.y, F2.u, F2.v, 'AutoScaleFactor', 20);
% quiver(F3.x, F3.y, F3.u, F3.v, 'AutoScaleFactor', 20);
% quiver(F4.x, F4.y, F4.u, F4.v, 'AutoScaleFactor', 20);
plot(segment5.x, segment5.y, 'r.-');
plot(segment6.x, segment6.y, 'r.-');
plot(segment7.x, segment7.y, 'r.-');
plot(segment8.x, segment8.y, 'r.-');
plot(segment1.x, segment1.y, 'r.-');
plot(segment2.x, segment2.y, 'r.-');
plot(segment3.x, segment3.y, 'r.-');
plot(segment4.x, segment4.y, 'r.-');
plot(segment9.x, segment9.y, 'b.-');
plot(segment10.x, segment10.y, 'b.-');
plot(segment11.x, segment11.y, 'b.-');
plot(segment12.x, segment12.y, 'b.-');
plot(segment13.x, segment13.y, 'g.-');
plot(segment14.x, segment14.y, 'g.-');
plot([Q1(1) Q1_(1)], [Q1(2) Q1_(2)], 'k--');
plot([Q2(1) Q2_(1)], [Q2(2) Q2_(2)], 'k--');
plot([Q3(1) Q3_(1)], [Q3(2) Q3_(2)], 'k--');
plot([Q4(1) Q4_(1)], [Q4(2) Q4_(2)], 'k--');
plot(V1(1), V1(2), 'g*');
plot(V2(1), V2(2), 'g*');
% plot(T1(1), T1(2), 'b*');
% plot(T2(1), T2(2), 'b*');
% plot(Q1(1), Q1(2), 'b*-');
% plot(segment5.x(1), segment5.y(1), 'b*-');
% plot(X4, Y4, 'g.-');

figure(2);
hold on;
grid on;
xlim([0, 900]);
ylim([0, 600]);
xlabel('x [m]');
ylabel('y [m]');
plot(s1.x, s1.y, 'r.-');
plot(s1.x(1), s1.y(1), 'b*');
plot(s2.x, s2.y, 'r.-');
plot(s2.x(1), s2.y(1), 'b*');
plot(s3.x, s3.y, 'r.-');
plot(s3.x(1), s3.y(1), 'b*');
plot(s4.x, s4.y, 'r.-');
plot(s4.x(1), s4.y(1), 'b*');
plot(s5.x, s5.y, 'r.-');
plot(s5.x(1), s5.y(1), 'b*');
plot(s6.x, s6.y, 'r.-');
plot(s6.x(1), s6.y(1), 'b*');
plot(s7.x, s7.y, 'r.-');
plot(s7.x(1), s7.y(1), 'b*');
plot(s8.x, s8.y, 'r.-');
plot(s8.x(1), s8.y(1), 'b*');
plot(s9.x, s9.y, 'r.-');
plot(s9.x(1), s9.y(1), 'b*');
plot(s10.x, s10.y, 'r.-');
plot(s10.x(1), s10.y(1), 'b*');
plot(s11.x, s11.y, 'r.-');
plot(s11.x(1), s11.y(1), 'b*');
plot(s12.x, s12.y, 'r.-');
plot(s12.x(1), s12.y(1), 'b*');
plot(s13.x, s13.y, 'r.-');
plot(s13.x(1), s13.y(1), 'b*');
plot(s14.x, s14.y, 'r.-');
plot(s14.x(1), s14.y(1), 'b*');

figure(3);
hold on;
grid on;
xlim([0, 900]);
ylim([0, 600]);
xlabel('x [m]');
ylabel('y [m]');
plot(map.x, map.y, 'b.-');
% quiver(335.8, 326.1,cos(phi4_0), sin(phi4_0), 'AutoScaleFactor', 20);

rmpath('functions');
rmpath('data');
