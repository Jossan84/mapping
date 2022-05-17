% testVehicleControl.m
% 05/06/2019
% Test vehicle control
% - Simulation Parameters
% - Set Models & Initialization
% - Make Map
% - Path Planning
% - High Level Control
% - Low  Level Control

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
addpath('path')

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
    n = 450; % [] (Number of simulation steps)
end

%% **************************** Set Models & Initialization **************************************
% *********************************** Initialization *********************************************
% Simulation parameters
T = 40e-3;% [s]
SetPoint_Vx = 20 / 3.6; %[m/s]
iBezier = 1;
tBezier = 0.0;
tOffset = 0.3; %To search goalPoint ahead to the car
maxDistance = 0.1;
maxIter = 25;

% Car Parameters
l             = 2.7;          % [m][Wheelbase]
steeringRatio = 15.8;       % []

% Path Plannig Parameters  
pointDistance = 10;         % [m]
LAMBDA = 1;
% Path Plannig Outputs
RefStates = struct('t', 0, 'x', 0, 'y', 0, 'v_x', 0, 'v_y', 0, 'yaw', 0, 'yawRate',0);

% ************************************ Set Models ************************************************
% Real car
dynamicCar   = DynamicModel('T',T,'m',1500,'I_z',1000,'C_f',20000,'C_r',30000,'l',l, ...
                            'Init',struct('t', 0, 'x', 5, 'y', 0, 'v_x', SetPoint_Vx,...
                            'v_y', 0, 'yaw', pi/2, 'yawRate',0));              
kinematicCar = KinematicModel('T',T,'l',l,'Init',struct('t', 0, 'x', 5, 'y', 0, 'v_x',...
                              SetPoint_Vx, 'v_y', 0, 'yaw', pi/2, 'yawRate',0));                              
realCar = dynamicCar;
%realCar  = kinematicCar;

% Actuation car
steeringActuator = ActuationModel('T',T,'steeringRatio',steeringRatio,'maxSteeringWheelAngle',400,...
                                  'minSteeringWheelAngle', -400);
% Virtual car
virtualKinematicCar = KinematicModel('T',T,'l',l,'Init',struct('t', 0, 'x', 5, 'y', 0, 'v_x',...
                                     SetPoint_Vx, 'v_y', 0, 'yaw', pi/2, 'yawRate',0));
virtualDynamicCar   = DynamicModel('T',T,'m',1500,'I_z',1000,'C_f',20000,'C_r',30000,'l',l, ...
                                   'Init',struct('t', 0, 'x', 5, 'y', 0, 'v_x', SetPoint_Vx,...
                                   'v_y', 0, 'yaw', pi/2, 'yawRate',0));
virtualCar  = virtualKinematicCar;

% Get Init Car States
realCarStates(1)    = realCar.getStates(realCar);
virtualCarStates(1) = virtualCar.getStates(virtualCar);

%% ************************************ LOAD PATH **************************************************
path_G = Path('Bezier3');

%% ********************************* SIMULATION LOOP **********************************************   
for k = 1 : n - 1
    %% ************************************ Get Path **********************************************
    
    % From Global to Local Path                                                      
    path_L = globalPath2localPath(path_G,realCar.getStates(realCar));
    
    % Get CTE; 
    [cte(k),iBezier,tBezier,nIter] = path_L.getCrossTrackError(path_L,0,0,0,iBezier,tBezier,...
                                                               maxDistance,maxIter,'Newton');
    
    %Get the Path that is Forward
    % TO DO
    
    %Check if the map ends
    % TO DO
%     if isempty(Path_L)
%       disp(['End of map: Exit k ' num2str(k)]);
%       plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates,...
%                      'Path_G', Path_G,'cte', cte, 'delta', [delta,0], 'T', T, 'iter', k-1,...
%                      'MarkerSize', markerSize,'RefStates',RefStates);
%       return;
%     end
    %% TEST
        if k == 304
            ;
        end
    %%
    % ********************************  Path Planning  *****************************************
    %Get the look ahead point
    if (test == 0)
        if(k == 1)                  
          [goalPoint, isFound, nIter1] = path_L.getLookAheadPoint(path_L,0,0,pointDistance,iBezier,...
                                                          tBezier+tOffset,maxDistance,maxIter,'Newton');                                                
          global firstGoalPoint
          firstGoalPoint = goalPoint - [realCarStates(1).y;realCarStates(1).x];
        end
    else        
          [goalPoint, isFound, nIter1] = path_L.getLookAheadPoint(path_L,0,0,pointDistance,iBezier,...
                                                          tBezier+tOffset,maxDistance,maxIter,'Newton');
          if(isFound ~= 1)
             disp('Not Found'); 
          end
    end   
         
    if (test == 0)
        if(k == 1)
            virtualCar  = virtualCar.setStates(virtualCar,realCar.getStates(realCar));
            deltaP      = atan2(2 * l * goalPoint(2), pointDistance^2);
            virtualCar  = virtualCar.update(virtualCar,deltaP,realCarStates(k).v_x); 
            RefStates   = virtualCar.getStates(virtualCar);
        end
    else
        virtualCar  = virtualCar.setStates(virtualCar,realCar.getStates(realCar));
            for i=1:LAMBDA
                
                % From Global to Local Path (Virtual Car)
                path_L = globalPath2localPath(path_G,virtualCar.getStates(virtualCar));
                
                % Get the Path that is Forward
                % TO DO
       
                % Check if the map ends
                % TO DO
%                 if isempty(Path_L)
%                     disp(['End of map: Exit k ' num2str(k)]);
%                     plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates,...
%                                    'Path_G', Path_G,'cte', cte, 'delta', [delta,0], 'T', T, 'iter', k-1,...
%                                    'MarkerSize', markerSize,'RefStates',RefStates);
%                     return;
%                 end
                                                                                     
                [goalPoint, isFound, nIter1] = path_L.getLookAheadPoint(path_L,0,0,pointDistance,iBezier,...
                                                          tBezier+tOffset,maxDistance,maxIter,'Newton');                                             
                                                            
                deltaP      = atan2(2 * l * goalPoint(2), pointDistance^2);
                virtualCar  = virtualCar.update(virtualCar,deltaP,realCarStates(k).v_x); 
                RefStates(i)= virtualCar.getStates(virtualCar);
            end
    end
          
    %% ******************************** Set point (Add Dinamic to PathPlanning)************************
    if false
%         nSetPoint = (goalPoint(1)^2 + goalPoint(2)^2) / (2 * goalPoint(2)) * atan2(2 * goalPoint(1) * goalPoint(2), goalPoint(1)^2 - goalPoint(2)^2) / (SetPoint_Vx * T);
        LAMBDA = 120;%floor(nSetPoint);

        a1 = 0.6;
    %    terma = 0;
    %    termb = 0;
    %    for i = 0 : LAMBDA - 1
    %        terma = terma + a1^i;
    %    end
    %    for i = 0 : LAMBDA - 2
    %        auxTerm = 0;
    %        for j = 0 : i
    %            auxTerm = auxTerm + a1^j;
    %        end
    %        termb = termb + auxTerm;
    %    end
        terma = (1 - a1^LAMBDA) / (1 - a1);
        termb = (LAMBDA - terma) / (1 - a1);
        
        v_y_b1 = (LAMBDA * RefStates.v_y - terma * currentState.v_y) / (RefStates.v_y * termb);
        yawRate_b1 = (LAMBDA * RefStates.yawRate - terma * currentState.yawRate) / (RefStates.yawRate * termb);
        
        v_y_d = zeros(LAMBDA + 1, 1);    
        yawRate_d = zeros(LAMBDA + 1, 1);
        
        v_y_d(1) = currentState.v_y;
        yawRate_d(1) = currentState.yawRate;
        
        for i = 2 : LAMBDA + 1
            v_y_d(i) = a1 * v_y_d(i - 1) + v_y_b1 * RefStates.v_y;
            yawRate_d(i) = a1 * yawRate_d(i - 1) + yawRate_b1 * RefStates.yawRate;
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
         Y   = vertcat(Y,([RefStates(1).v_y; RefStates(1).yawRate]...
                         - ((system.Ad^j) * [currentState.v_y;currentState.yawRate])));
%        Y   = vertcat(Y,([v_y_d(j + 1); yawRate_d(j + 1)]...
%                         - ((system.Ad^j) * [currentState.v_y;currentState.yawRate])));
%        Y   = vertcat(Y,([v_y_d(j + 1); RefStates.yawRate]...
%                         - ((system.Ad^j) * [currentState.v_y;currentState.yawRate])));
        sum = sum + system.Ad^(j-1);
        X   = vertcat(X,(sum * system.Bd));
    end
    
    % Compute linear regression to obtain the delta that minimizes the path trajectory and
    % get the steering wheel angle. 
    u(k) = Y'/X';
    
    steeringWheelAngle = u(k) * steeringRatio; % Desired steering wheel angle [rad]    

    
    %% ************************** Car/Actuation [EPS] ********************************
    
    v_x(k)            = SetPoint_Vx;  % Desired longitudinal speed [m/s]
    %steeringActuator  = steeringActuator.update(steeringActuator,steeringWheelAngle);
    %delta(k)          = steeringActuator.states.delta;
    delta(k) = steeringWheelAngle / steeringRatio;
        
    % Update car step movement
    realCar = realCar.update(realCar,delta(k),v_x(k));
    
    %Get states to plot
    realCarStates(k+1)    = realCar.getStates(realCar);
    virtualCarStates(k+1) = virtualCar.getStates(virtualCar);
    %desiredStates(:,1)    = [0,0];
    %desiredStates(:,k+1)  = [RefStates.v_y; RefStates.yawRate];
    
end
%% ************************************ END LOOP *************************************************

%% ************************************** Plots **************************************************
plotSimulation('realCarStates', realCarStates, 'virtualCarStates', virtualCarStates, 'Path_G',...
                path_G,'cte', [cte 0], 'delta', [delta, 0], 'T', T, 'iter', k,...
                'MarkerSize', markerSize,'RefStates',RefStates);

error = cte * cte';
disp(['Error: ', num2str(error)]);

rmpath('functions');
rmpath('models');
rmpath('map');
rmpath('path');
