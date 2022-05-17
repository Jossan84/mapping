%getCrossTrackError
%06/05/2019

%% function [cte] = getCrossTrackError(Path_L)
% Description:
% This function return the cross track error between the car (0,0) and the path.
% The method for this calculation is explained in the next link:
% [https://es.wikipedia.org/wiki/Distancia_de_un_punto_a_una_recta]

% Input:
%   Path_L = Local Path [m]
% Output:
%   cte    = Cross Track Error [m]

function [cte] = getCrossTrackError(Path_L)
    
    [n m] = size(Path_L);
    
    %Get the index of the nearest point
    [~, minIndex] = min(pointDistances(Path_L,0).^2);
    
    if minIndex == m
        cte = NaN;
        return;
    else
        %Choose the points for search min cte
        % p1 --> Nearest point
        % p2 --> Next point to p1
        % p2 --> Previous point to p1
    
        p1.x = Path_L(1,minIndex);
        p1.y = Path_L(2,minIndex);
        if minIndex-1 == 0
            p0 = p1;
        else
            p0.x = Path_L(1,minIndex-1);
            p0.y = Path_L(2,minIndex-1);
        end
        p2.x = Path_L(1,minIndex+1);
        p2.y = Path_L(2,minIndex+1);    
    
        %Get Cross Track Error
        %r1    --> Projection of car point (0,0) on the straight from p0 to p1
        %r2    --> Projection of car point (0,0) on the straight from p1 to p2
        %norm1 --> Morm of p0 and p1
        %norm2 --> Morm of p1 and p2
        %d1    --> Min distance from car to the straight p0 to p1
        %d2    --> Min distance from car to the straight p1 to p2
    
        r1    = ((p1.y*0  - p0.y*0 - p1.y*p0.x + p0.y*p0.x) -  ...
                ( p1.x*0  - p0.x*0 - p1.x*p0.y + p0.x*p0.y));
        r2    = ((p2.y*0  - p1.y*0 - p2.y*p1.x + p1.y*p1.x) -  ...
                ( p2.x*0  - p1.x*0 - p2.x*p1.y + p1.x*p1.y));    
        norm1 = sqrt((p1.x-p0.x)^2 + (p1.y-p0.y)^2);
        norm2 = sqrt((p2.x-p1.x)^2 + (p2.y-p1.y)^2);
    
    
        %Prevent division by zero or by NaN
        if norm1 == 0 && r1 == 0
        d1 = NaN; 
        else    
        d1 = r1/norm1;   
        end
    
        if norm2 == 0 && r2 == 0
        d2 = NaN; 
        else    
        d2 = r2/norm2;   
        end
       
        cte = min([d1, d2]);
    end
    
end

