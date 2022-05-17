function actuationModel = ActuationModel(varargin)
    
    % Default params
    params.T                       = 40e-3;         % [s]
    params.steeringRatio           = 15.8 ;         % []
    params.maxSteeringWheelAngle   = deg2rad(600);  % [rad]
    params.minSteeringWheelAngle   = deg2rad(-600); % [rad]
        
   % Defaul states
%     persistent y_k;
%     y_k = 0;    
    states = struct('t', 0, 'delta',0,'y_k',0);
    
    for i = 1 :2: nargin
        switch varargin{i}
            case 'T'
                params.T = varargin{i + 1};           
            case 'steeringRatio'
                params.steeringRatio = varargin{i + 1}; 
            case 'maxSteeringWheelAngle'
                params.maxSteeringWheelAngle = deg2rad(varargin{i + 1});
            case 'minSteeringWheelAngle'
                params.minSteeringWheelAngle = deg2rad(varargin{i + 1});  
            otherwise
                error('Wrong argument');
        end
    end
            

    actuationModel.params = params;
    actuationModel.states = states;
    
    actuationModel.update    = @update;
    actuationModel.getStates = @getStates;
    actuationModel.setStates = @setStates;
end

function [actuationModel] = update(actuationModel,steeringWheelAngle)
    
    maxDelta = actuationModel.params.maxSteeringWheelAngle / actuationModel.params.steeringRatio;
    minDelta = actuationModel.params.minSteeringWheelAngle / actuationModel.params.steeringRatio;
    
    delta = steeringWheelAngle / actuationModel.params.steeringRatio;
    
    % Actuation dynamic
    a = 0.75;
    b = 0.25;

    actuationModel.states.y_k = a*actuationModel.states.y_k + b*steeringWheelAngle;
    
    delta = actuationModel.states.y_k / actuationModel.params.steeringRatio;
    
    % Saturation for max a min angle of delta
    if (delta >= maxDelta)
        actuationModel.states.delta = maxDelta;
    elseif (delta < minDelta)
        actuationModel.states.delta = minDelta;
    else
        actuationModel.states.delta = delta;
    end
    
    actuationModel.states.t   = actuationModel.states.t + actuationModel.params.T;
    
end

function states = getStates(actuationModel)
    states = actuationModel.states;
end

function actuationModel = setStates(actuationModel,t,delta,y_k)
    actuationModel.states = struct('t', t, 'delta',delta,'y_k',y_k);
end
