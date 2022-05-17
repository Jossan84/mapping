% testVehicleControl.m
% 05/06/2019
% Test vehicle control
% - Simulation Parameters
% - Set Models & Initialization
% - Navigation
% - Local Path Planning
% - Control
% - Actuation


close all;
clear;
more off;
clc;

if exist('OCTAVE_VERSION', 'builtin') ~= 0% OCTAVE
    IDE = 'OCTAVE';
    markerSize = 12;
else% MATLAB
    IDE = 'MATLAB';
    markerSize = 6;
end

addpath('functions');
addpath('models');
addpath('map');

prompt = {'Options: [0]Go one Point [1]Path Following'};
dlgtitle = 'Choose option';
definput = {'0'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,[1 60],definput,opts);

global test;
if (str2num(answer{1}) == 0)
    test = 0;
    n = 200; % [] (Number of simulation steps)
else
    test = 1;
    n = 3000; % [] (Number of simulation steps)
end

%% **************************** Set Models & Initialization **************************************
% *********************************** Initialization *********************************************
% Simulation parameters
T = 40e-3;% [s]
SetPoint_Vx = 30 / 3.6; %[m/s]

% Initial values
% xInit = 5;
% yInit = 0;
% yawInit = pi/2;
xInit = 444.6;
yInit = 234.8;
yawInit = deg2rad(320);

% Car Parameters
l             = 2.7;          % [m][Wheelbase]
steeringRatio = 15.8;       % []

% Path Plannig Parameters  
pointDistance = 10;         % [m]
LAMBDA = 1;
A = [0,0;0,0.6];
B = [0,0.4];
% Path Plannig Outputs
RefStates = struct('t', 0, 'x', 0, 'y', 0, 'v_x', 0, 'v_y', 0, 'yaw', 0, 'yawRate',0);

% ************************************ Set Models ************************************************
% Real car
dynamicCar   = DynamicModel('T',T,'m',1500,'I_z',1000,'C_f',20000,'C_r',30000,'l',l, ...
                            'Init',struct('t', 0, 'x', xInit, 'y', yInit, 'v_x', SetPoint_Vx,...
                            'v_y', 0, 'yaw', yawInit, 'yawRate',0));              
kinematicCar = KinematicModel('T',T,'l',l,'Init',struct('t', 0, 'x', xInit, 'y', yInit, 'v_x',...
                              SetPoint_Vx, 'v_y', 0, 'yaw', yawInit, 'yawRate',0));
GeometricCar = GeometricModel('T',T,'l',l,'Init',struct('t', 0, 'x', xInit, 'y', yInit, 'v_x',...
                                     SetPoint_Vx, 'v_y', 0, 'yaw', yawInit, 'yawRate',0));
realCar = dynamicCar;
%realCar  = kinematicCar;


% Actuation car
steeringActuator = ActuationModel('T',T,'steeringRatio',steeringRatio,'maxSteeringWheelAngle',400,...
                                  'minSteeringWheelAngle', -400);
% Virtual car
virtualKinematicCar = KinematicModel('T',T,'l',l,'Init',struct('t', 0, 'x', xInit, 'y', yInit, 'v_x',...
                                     SetPoint_Vx, 'v_y', 0, 'yaw', yawInit, 'yawRate',0));
virtualGeometricCar = GeometricModel('T',T,'l',l,'Init',struct('t', 0, 'x', xInit, 'y', yInit, 'v_x',...
                                     SetPoint_Vx, 'v_y', 0, 'yaw', yawInit, 'yawRate',0));
virtualDynamicCar   = DynamicModel('T',T,'m',1500,'I_z',1000,'C_f',20000,'C_r',30000,'l',l, ...
                                   'Init',struct('t', 0, 'x', xInit, 'y', yInit, 'v_x', SetPoint_Vx,...
                                   'v_y', 0, 'yaw', yawInit, 'yawRate',0));
%virtualCar  = virtualKinematicCar;                               
virtualCar  = virtualGeometricCar;


% Get Init Car States
realCarStates(1)    = realCar.getStates(realCar);
virtualCarStates(1) = virtualCar.getStates(virtualCar);

%% ************************************ LOAD MAP **************************************************
% load('Map.mat');
load('loopMap.mat');
Path_G.x    = Map.x;
Path_G.y    = Map.y;

%% ********************************* SIMULATION LOOP **********************************************   
for k = 1 : n - 1
    %% ******************************* Navigation *************************************************
    % TO DO: Select the route and load the map segments for this route.
    % Notes: Here we choose te destination point and obtain the segments of
    % the map to achieve this point. If we need to recalculate the route
    % here is the place where is it done and if we have to park at the
    % destination here we take into account the segments that lead us to
    % the park slot (For this we have to include the parking map into the maps).
    % Note: We can reduce the error of the map by the use of the local
    % sensor of the car (camera, liddar, radar, etc...).
    
    %% ************************************ Get Local Path ****************************************
    % From Global to Local (Only for ger the CTE);
    [Path_L]= global2local(Path_G, realCar.getStates(realCar));
    cte(k)  = getCrossTrackError(Path_L);
           
    %Get the Path that is Forward
    Path_L  = Path_L(:, Path_L(1, :) > -1e-10);
       
    %Check if the map ends
    if isempty(Path_L)
      disp(['End of map: Exit k ' num2str(k)]);
      plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates,...
                     'Path_G', Path_G,'cte', cte, 'delta', [delta,0], 'T', T, 'iter', k-1,...
                     'MarkerSize', markerSize,'RefStates',RefStates);
      return;
    end
    
    %% ********************************  Local Path Planning  *****************************************
    %Get the look ahead point
    if (test == 0)
        if(k == 1)
          [goalPoint  ,minDistIndex]  = getLookAheadPoint('IDE', IDE, 'Path_L', Path_L,...
                                                          'pointDistance', pointDistance);
          global firstGoalPoint
          firstGoalPoint = goalPoint - [realCarStates(1).y;realCarStates(1).x];
        end
    else
          [goalPoint  ,minDistIndex]  = getLookAheadPoint('IDE', IDE, 'Path_L', Path_L,...
                                                          'pointDistance', pointDistance);    
    end   
         
    if (test == 0)
        if(k == 1)
%             virtualCar  = virtualCar.setStates(virtualCar,realCar.getStates(realCar));
%             deltaP      = atan2(2 * l * goalPoint(2), pointDistance^2);
%             virtualCar  = virtualCar.update(virtualCar,deltaP,realCarStates(k).v_x); 
%             RefStates   = virtualCar.getStates(virtualCar);
            virtualCar  = virtualCar.setStates(virtualCar,realCar.getStates(realCar));
            for i=1:LAMBDA
                [Path_L]= global2local(Path_G, virtualCar.getStates(virtualCar));
                % Get the Path that is Forward
                Path_L  = Path_L(:, Path_L(1, :) > -1e-10);
       
                % Check if the map ends
                if isempty(Path_L)
                    disp(['End of map: Exit k ' num2str(k)]);
                    plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates,...
                                   'Path_G', Path_G,'cte', cte, 'delta', [delta,0], 'T', T, 'iter', k-1,...
                                   'MarkerSize', markerSize,'RefStates',RefStates);
                    return;
                end
                
                [goalPoint  ,minDistIndex]  = getLookAheadPoint('IDE', IDE, 'Path_L', Path_L,...
                                                                'pointDistance', pointDistance);

                deltaP       = atan2(2 * l * goalPoint(2), pointDistance^2);
   
                yawRateRef(i)   = ((2 * l * goalPoint(2))/pointDistance^2) * (realCarStates(k).v_x/l);
                
                % Add Dynamic to yawRate
                if i == 1
                    yawRateDes(i)  = firstOrderTransferFunction(yawRateRef(i),0,A(2,2),B(2)); 
                else
                    yawRateDes(i)  = firstOrderTransferFunction(yawRateRef(i),yawRateDes(i-1),A(2,2),B(2));
                end
                                                
                % Distance to Reference
                [s,t,n] = time2Reference(goalPoint,pointDistance,realCarStates(k).v_x,T);
                                
                % Gain
                Y = firstOrderDynamicsOutputSum(yawRateDes(i),yawRateRef(i),n,A(2,2),B(2));
                K(i) = (n+1)*(yawRateRef(i)- yawRateDes(i))/Y;
                
                DesStates(i) = virtualCar.getStates(virtualCar);
                DesStates(i).yawRate = K(i) * yawRateDes(i);
                virtualCar  = virtualCar.setStates(virtualCar,DesStates(i));
                
                virtualCar   = virtualCar.update(virtualCar,deltaP,realCarStates(k).v_x);
                RefStates(i) = virtualCar.getStates(virtualCar);
                                
            end
        end
    else
        virtualCar  = virtualCar.setStates(virtualCar,realCar.getStates(realCar));
            for i=1:LAMBDA
                [Path_L]= global2local(Path_G, virtualCar.getStates(virtualCar));
                % Get the Path that is Forward
                Path_L  = Path_L(:, Path_L(1, :) > -1e-10);
       
                % Check if the map ends
                if isempty(Path_L)
                    disp(['End of map: Exit k ' num2str(k)]);
                    plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates,...
                                   'Path_G', Path_G,'cte', cte, 'delta', [delta,0], 'T', T, 'iter', k-1,...
                                   'MarkerSize', markerSize,'RefStates',RefStates);
                    return;
                end
                
                [goalPoint  ,minDistIndex]  = getLookAheadPoint('IDE', IDE, 'Path_L', Path_L,...
                                                                'pointDistance', pointDistance);
                                                            
                deltaP       = atan2(2 * l * goalPoint(2), pointDistance^2);
   
                yawRateRef(i)   = ((2 * l * goalPoint(2))/pointDistance^2) * (realCarStates(k).v_x/l);
                
                % Add Dynamic to yawRate
                if i == 1
                    yawRateDes(i)  = firstOrderTransferFunction(yawRateRef(i),0,A(2,2),B(2)); 
                else
                    yawRateDes(i)  = firstOrderTransferFunction(yawRateRef(i),yawRateDes(i-1),A(2,2),B(2));
                end
                                                
                % Distance to Reference
                [s,t,n] = time2Reference(goalPoint,pointDistance,realCarStates(k).v_x,T);
                                
                % Gain
                Y = firstOrderDynamicsOutputSum(yawRateDes(i),yawRateRef(i),n,A(2,2),B(2));
                K(i) = (n+1)*(yawRateRef(i)- yawRateDes(i))/Y;
                
                DesStates(i) = virtualCar.getStates(virtualCar);
                DesStates(i).yawRate = K(i) * yawRateDes(i);
                virtualCar  = virtualCar.setStates(virtualCar,DesStates(i));
                
                virtualCar   = virtualCar.update(virtualCar,deltaP,realCarStates(k).v_x);
                RefStates(i) = virtualCar.getStates(virtualCar);
                                
            end
    end
          
    %% ******************************** Control *****************************************
    %Initialization
    Y   = [];
    X   = [];
    sum = zeros(2);
    %Get system 
    system = realCar.getSystem(realCar);
    currentState  = realCar.getStates(realCar);

    for j=1:LAMBDA
%        Y   = vertcat(Y,([RefStates(j).v_y; RefStates(j).yawRate]...
%                          - ((system.Ad^j) * [currentState.v_y;currentState.yawRate])));
%        Y   = vertcat(Y,([DesStates(j).v_y; DesStates(j).yawRate]...
%                         - ((system.Ad^j) * [currentState.v_y;currentState.yawRate])));
       Y   = vertcat(Y,([0; DesStates(j).yawRate]...
                        - ((system.Ad^j) * [0;currentState.yawRate])));               
       sum = sum + system.Ad^(j-1);
       X   = vertcat(X,(sum * system.Bd));
        
    end

    % Compute linear regression to obtain the delta that minimizes the path trajectory and
    % get the steering wheel angle. 
    u(k) = Y'/X';
    %u(k) = deltaP;
    steeringWheelAngle = u(k) * steeringRatio; % Desired steering wheel angle [rad]    

    
    %% ************************** Car/Actuation [EPS] ********************************
    
    v_x(k)            = SetPoint_Vx;  % Desired longitudinal speed [m/s]
%     steeringActuator  = steeringActuator.update(steeringActuator,steeringWheelAngle);
%     delta(k)          = steeringActuator.states.delta;
    delta(k) = steeringWheelAngle / steeringRatio;
        
    % Update car step movement
    realCar = realCar.update(realCar,delta(k),v_x(k));
    
    % Get states to plot
    realCarStates(k+1)    = realCar.getStates(realCar);
    virtualCarStates(k+1) = virtualCar.getStates(virtualCar);
    % Plot movement
%     figure(1);
%     ax = gca;
%     plotStep(Map,realCarStates,k,ax);
    
end

steeringActuator = steeringActuator.setStates(steeringActuator,0,0,0);
%% ************************************ END LOOP *************************************************

%% ************************************** Plots **************************************************
plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates, 'Path_G',...
                Path_G,'cte', [cte 0], 'delta', [delta, 0], 'T', T, 'iter', k,...
                'MarkerSize', markerSize,'RefStates',RefStates);

error = cte * cte';
disp(['Error^2: ', num2str(error)]);

rmpath('functions');
rmpath('models');
rmpath('map');
