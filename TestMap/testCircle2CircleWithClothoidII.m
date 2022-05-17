%testCircle2CircleWithClothoidII
%23/11/2020

close all;
clear;
clc

addpath('functions');

%% Description
% Connection between a oriented circle with another oriented circle.

% Circles center points
C1 = [-3, 3];
F1 = [-4, 3];
% C1 = [-3, 3.61];
% F1 = [-4, 3.61];
R1 = sqrt((F1(1)-C1(1))^2 + ((F1(2)-C1(2))^2));
lambdaC1 = 1;
yaw1 = linspace(0, 2*pi, 100);
x1 = C1(1) + (R1 * cos(yaw1));
y1 = C1(2) + (R1 * sin(yaw1));

C2 = [1, 1];
F2 = [3, 1];
R2 = sqrt((F2(1)-C2(1))^2 + ((F2(2)-C2(2))^2));
lambdaC2 = -1;
yaw2 = linspace(0, 2*pi, 100);
x2 = C2(1) + (R2 * cos(yaw2));
y2 = C2(2) + (R2 * sin(yaw2));

% Directions
phi1_0 = pi * 7/8;
phi2_0 = phi1_0 + pi;

[X1, Y1] = getClothoidPoints(F1, C1, lambdaC1, phi1_0, pi * 3/2);
[X2, Y2] = getClothoidPoints(F2, C2, lambdaC2, phi2_0, pi * 3/2);

% Points
P1 = [X1(1), Y1(1)];
 
startPoint_Y = P1(2);
startPoint_X = P1(1);
endPoint_X = 6;
m = tan(phi1_0);
n =  startPoint_Y - (tan (phi1_0) * startPoint_X );
endPoint_Y = m*endPoint_X + n;
P1_ = [endPoint_X, endPoint_Y];

P2 = P1_;
startPoint_Y = P2(2);
startPoint_X = P2(1);
endPoint_X = -6;
m = tan(phi2_0);
n =  startPoint_Y - (tan (phi2_0) * startPoint_X );
endPoint_Y = m*endPoint_X + n;
P2_ = [endPoint_X, endPoint_Y];


line1.x1 = P1(1);
line1.y1 = P1(2);
line1.x2 = P1_(1);
line1.y2 = P1_(2);
line2.x1 = P2(1);
line2.y1 = P2(2);
line2.x2 = P2_(1);
line2.y2 = P2_(2);

V = getIntersectionPoint(line2, line1);
phi1_0 = deg2rad(getAngle(P1, V));
phi2_0 = deg2rad(getAngle(P2, V));

% Plot
figure('Name', 'Plot Clothoid');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
% xlim(ax, [-2 8]);
ylim(ax, [-2 8]);
xlabel(ax, 'x [m]');
ylabel(ax, 'y [m]');
plot(ax, x1, y1, 'b.-');
plot(ax, x2, y2, 'b.-');
plot(ax, F1(1), F1(2), 'ro');
plot(ax, F2(1), F2(2), 'ro');
plot(ax, X1, Y1, 'r.-');
plot(ax, X2, Y2, 'r.-');
% plot(ax, X3, Y3, 'r.-');
% plot(ax, X4, Y4, 'r.-');
% plot(ax, P1_(1), P1_(2), 'k*');
% plot(ax, P2_(1), P2_(2), 'k*');
plot(ax, [P1(1) P1_(1)], [P1(2) P1_(2)], 'b--');
plot(ax, [P2(1) P2_(1)], [P2(2) P2_(2)], 'b--');
% plot(ax, [P1(1), V(1)], [P1(2), V(2)], 'b--');
% plot(ax, [P2(1), V(1)], [P2(2), V(2)], 'b--');
plot(ax, P1(1), P1(2), 'k*');
plot(ax, P1_(1), P1_(2), 'k*');
plot(ax, P2(1), P2(2), 'k*');
plot(ax, P2_(1), P2_(2), 'k*');
% plot(ax, T1(1), T1(2), 'b*');
% plot(ax, T2(1), T2(2), 'b*');
text(ax, P1(1), P1(2)+0.2, 'P1','FontSize',14);
text(ax, P2(1), P2(2)+0.2, 'P2','FontSize',14);
text(ax, V(1), V(2) + 0.2, 'V','FontSize',14);
plot(ax, V(1), V(2), 'g*');
quiver(ax, F1(1), F1(2), cos(pi * 3/2), sin(pi * 3/2), 'Color','k', 'AutoScaleFactor', 1);
quiver(ax, F2(1), F2(2), cos(pi * 3/2), sin(pi * 3/2), 'Color','k', 'AutoScaleFactor', 1);

rmpath('functions');