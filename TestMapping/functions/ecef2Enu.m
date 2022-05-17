function [xEast,yNorth,zUp] = ecef2Enu(X, Y, Z, lat0, lon0)
% ecef2Enu
% 10/11/2021
% Transform geocentric Earth-centered Earth-fixed coordinates to local east-north-up
%
% Usage:
%     [xEast,yNorth,zUp] = ecef2Enu(X, Y, Z, lat0, lon0, lat, lon)
%     arguments:
%         lat0      	 Reference geodetic latitude [rad]
%         lon0           Reference geodetic longitude [rad]
%         X              ECEF x coordinates [m]
%         Y              ECEF y coordinates [m]
%         Z              ECEF z coordinates [m]
%
% Examples:
%     [xEast, yNorth, zUp] = ecef2Enu(X, Y, Z, lat0, lon0)
%
% Dependencies:
%     
  %% Vectorial method
  xEast = -sin(lon0) .* (X - X(1)) + cos(lon0) .* (Y - Y(1));
  yNorth = -sin(lat0).*cos(lon0) .* (X - X(1)) + -sin(lat0).*sin(lon0) .* (Y - Y(1)) + cos(lat0) .* (Z - Z(1));
  zUp = cos(lat0).*cos(lon0) .* (X - X(1)) + cos(lat0).*sin(lon0) .* (Y - Y(1)) + sin(lat0) .* (Z - Z(1));
 
%% Iterative method  
%   R = [         -sin(lon0),            cos(lon0),          0;
%        -cos(lon0)*sin(lat0), -sin(lon0)*sin(lat0),  cos(lat0);
%         cos(lon0)*cos(lat0),  sin(lon0)*cos(lat0), sin(lat0)];
%     
%   xEast = zeros(size(X));
%   yNorth = zeros(size(Y));
%   zUp = zeros(size(Z));
%   for i=1:length(X)  
%       p = R * [X(i) - X(1); Y(i) - Y(1); Z(i) - Z(1)]; 
%       xEast(i,1) = p(1);
%       yNorth(i,1) = p(2);
%       zUp(i,1) = p(3);
%   end
%%  
end