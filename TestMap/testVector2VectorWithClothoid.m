%testVector2VectorWithClothoid
%16/11/2020

close all;
clear;
clc

addpath('functions');

%% Description
% Connection between a vector with direction and sense with another
% vector.

% Points
P1 = [0, 0];
P2 = [0, 10];
% V  = [8, 5];
V  = [2, 5];

P1V = V - P1;
P2V = V - P2;

phi1_0 = deg2rad(getAngle(P1, V));
phi2_0 = 2*pi - deg2rad(getAngle(P2, V));
beta = acos( (P1V(1)*P2V(1) + P1V(2)*P2V(2)) / (sqrt(P1V(1)^2 + P1V(2)^2) * sqrt(P2V(1)^2 + P2V(2)^2)));
lambda1 = 1;
lambda2 = -1;
omega = deg2rad(0);
tau = (pi - (omega + beta))/2;
R = 2;
L = 2*R*tau;
A = sqrt(R*L);
N = 100;
s = linspace(0, L, N);
ds = L/N;

termXL = 0;
termYL = 0;
for k = 1 : N
    termXL = termXL + cos(s(k)^2 / (2*A^2));
    termYL = termYL + sin(s(k)^2 / (2*A^2));
end

XL = ds * termXL;
YL = ds * termYL;

distT1J1 = XL;
distJ1H1 = YL * tan(tau);
distJ2H2 = distJ1H1;
distH1V = (R + YL/cos(tau)) * (sin(omega/2)/sin(beta/2));
distH2V = distH1V;
distT1V = distT1J1 + distJ1H1 + distH1V;
distT2V = distT1V;

T1 = V - distT1V * P1V/norm(P1V);
T2 = V - distT2V * P2V/norm(P2V);

x1 = zeros(length(s),1);
y1 = zeros(length(s),1);
x1(1) = T1(1);
y1(1) = T1(2);
x2 = zeros(length(s),1);
y2 = zeros(length(s),1);
x2(1) = T2(1);
y2(1) = T2(2);


for i = 1: length(s)
    x1(i+1) = x1(i) + ds * cos((lambda1*(s(i)^2)/(2*A^2)) + phi1_0);
    y1(i+1) = y1(i) + ds * sin((lambda1*(s(i)^2)/(2*A^2)) + phi1_0);
    
    x2(i+1) = x2(i) + ds * cos((lambda2*(s(i)^2)/(2*A^2)) + phi2_0);
    y2(i+1) = y2(i) + ds * sin((lambda2*(s(i)^2)/(2*A^2)) + phi2_0);
end

% Plot
figure('Name', 'Plot Clothoid');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
% xlim(ax, [-2 8]);
% ylim(ax, [-2 8]);
xlabel(ax, 'x [m]');
ylabel(ax, 'y [m]');
plot(ax, x1, y1, 'r.-');
plot(ax, x2, y2, 'r.-');
plot(ax, [P1(1), V(1)], [P1(2), V(2)], 'b--');
plot(ax, [P2(1), V(1)], [P2(2), V(2)], 'b--');
plot(ax, P1(1), P1(2), 'k*');
plot(ax, P2(1), P2(2), 'k*');
plot(ax, T1(1), T1(2), 'b*');
plot(ax, T2(1), T2(2), 'b*');
text(ax, P1(1), P1(2)+0.2, 'P1','FontSize',14);
text(ax, P2(1), P2(2)+0.2, 'P2','FontSize',14);
text(ax, V(1)+0.2, V(2), 'V','FontSize',14);
plot(ax, V(1), V(2), 'r*');

rmpath('functions');
