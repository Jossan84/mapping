%testVector2CircleWithClothoid
%16/11/2020

close all;
clear;
clc

addpath('functions');

%% Description
% Connection between a vector with direction and sense with a
% point of a circle.

% Vector
phi0 = 15*pi/8;
v = [cos(phi0); sin(phi0)];
% Circle
F = [3; 1];
C = [1; 1];
alpha = deg2rad(270);
lambdaC = -1; %[-1] clockwise [1] counterclockwise

% Get clothoid points
[x, y] = getClothoidPoints(F, C, lambdaC, phi0, alpha);

% Data for plotting
vT = lambdaC * [C(2)-F(2); C(1)-F(1)];
R = sqrt((F(1)-C(1))^2 + ((F(2)-C(2))^2));
[xC, yC] = circle(C(1), C(2), R);

% Plot
figure('Name', 'Plot Clothoid');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
xlim(ax, [-2 8])
ylim(ax, [-2 5])
xlabel(ax, 'x [m]');
ylabel(ax, 'y [m]');
plot(ax, xC, yC, 'b.-');
plot(ax, x, y, 'r.-');
text(ax, F(1)+0.1, F(2), 'F','FontSize',14);
text(ax, C(1), C(2)+0.2, 'C','FontSize',14);
plot(ax, [C(1),F(1)], [C(2),F(2)], 'b--');
quiver(ax, F(1),F(2),vT(1),-vT(2), 'Color','k', 'AutoScaleFactor',1);
quiver(ax, x(1),y(1),v(1),v(2), 'Color','k', 'AutoScaleFactor',6);


rmpath('functions');
