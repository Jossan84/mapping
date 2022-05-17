%plotSimulation
%05/08/2019

function plotSimulation( varargin )

    for i = 1 :2: nargin
        switch varargin{i}
            case 'realCarStates'
                realCarStates = varargin{i + 1};
            case 'virtualCarStates'
                virtualCarStates = varargin{i + 1};
            case 'Path_G'
                Path_G = varargin{i + 1};
            case 'cte'
                cte = varargin{i + 1};
            case 'delta'
                delta = varargin{i + 1};
            case 'T'
                T = varargin{i + 1};
            case 'iter'
                iter = varargin{i + 1};
            case 'MarkerSize'
                markerSize = varargin{i + 1};
            case 'RefStates'
                RefStates = varargin{i + 1}; 
            otherwise
                error('Wrong argument');
        end
    end
    
index = (1 : iter + 1)';
t = T * (index - 1);
global test
global firstGoalPoint
%% Plot inputs
figure(1);
plot([realCarStates.t], [realCarStates.v_x], 'b.-', 'MarkerSize', markerSize);
title('Longitudinal velocity');
ylabel('v_x [m/s]');
xlabel('Time [s]');
grid on;

figure(2);
plot(t, rad2deg(delta), 'b.-', 'MarkerSize', markerSize);
title('Wheel angle');
ylabel('delta [deg]');
xlabel('Time [s]');
grid on;

%% Plot outputs
figure(3);
plot([realCarStates.t], [realCarStates.v_y], 'b.-', 'MarkerSize', markerSize);
if(test == 0)
    hold on;
    plot(t,(RefStates.v_y),'r.-');
end
%hold on;
%plot(t,[RefStates.v_y],'r.-');
title('Lateral velocity');
ylabel('v_y [m/s]');
xlabel('Time [s]');
grid on;

figure(4);
plot([realCarStates.t], [realCarStates.yawRate], 'b.-', 'MarkerSize', markerSize);
if(test == 0)
    hold on;
    plot(t,(RefStates.yawRate),'r.-');
end
%hold on;
%plot(t,[RefStates.yawRate],'r.-')
title('Yaw rate');
ylabel('yaw rate [rad/s]');
xlabel('Time [s]');
grid on;

figure(5);
plot(t, cte, 'b.-', 'MarkerSize', markerSize);
title('Cross Track Error');
ylabel('cte [m]');
xlabel('Time [s]');
grid on;

figure(6);
plot([realCarStates.t], rad2deg([realCarStates.v_y]./[realCarStates.v_x]), 'b.-', 'MarkerSize', markerSize);
title('Slip Angle');
ylabel('Beta [deg]');
xlabel('Time [s]');
grid on;

% Plot xy
figure;
ax = gca;
axis(ax, 'equal');
hold(ax, 'on');
plot(ax, [realCarStates.x], [realCarStates.y], 'b.-', 'MarkerSize', markerSize);
hold on;
%plot(ax,[virtualCarStates.x], [virtualCarStates.y], 'r.-', 'MarkerSize', markerSize);
Path_G.plotPath(Path_G,1000);
% plot(Path_G.x,Path_G.y, 'k.-', 'MarkerSize', markerSize);
if(test == 0)
    plot(firstGoalPoint(2),firstGoalPoint(1), 'g*', 'MarkerSize', markerSize);
    ylim([0 100]);
    xlim([-20 20]);
end
xlabel('x [m]');
ylabel('y [m]');

% Plot Local Path
% figure(5);
% for i=1:n-1
% 
% plot(-plot_Goal_point(i,2),plot_Goal_point(i,1),'b*');
% hold on;
% plot(0,0,'r*');
% axis([-5 5 -5 15]);
% pause(0.8);
% hold off;

% figure(6);
% for i=1:n-1
%     hold on
%     for j=1:LAMBDA        
%         plot( RefStates(i,j).x,RefStates(i,j).y,'b*')
%     end
%     pause(0.1);
%     clf    
% end

end


