function [ Path_L] = global2local( Path_G, Position )

    % Rotation Matrix
    R_GL = [cos(Position.yaw) , sin(Position.yaw);
            -sin(Position.yaw), cos(Position.yaw)];
        
    d    = [Position.x; Position.y];  

    r    = -R_GL*d;                       % Tranlation
    Rot  = R_GL * [Path_G.x; Path_G.y];   % Rotation
    
    [n m] = size(Rot);
    Trans(1,:) = r(1) .* ones(1,m);
    Trans(2,:) = r(2) .* ones(1,m);
    
    Path_L    = Trans + Rot;    % Transformation Homogeneous
    
end
