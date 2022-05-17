function [X, Y, Z] = enu2Ecef(xEast, yNorth, zUp, lat0, lon0, h0)
% enu2Ecef
% 11/11/2021
% Transform geocentric Earth-centered Earth-fixed coordinates to local east-north-up
%
% Usage:
%     [X, Y, Z] = ecef2Enu(xEast, yNorth, zUp, lat0, lon0, h0)
%     arguments:
%         lat0      	   Reference geodetic latitude [rad]
%         lon0           Reference geodetic longitude [rad]
%         xEast          x ENU coordinates [m]
%         yNorth         y ENU coordinates [m]
%         zUp            z ENU coordinates [m]
%
% Examples:
%     [X, Y, Z] = ecef2Enu(xEast, yNorth, zUp, lat0, lon0, h0)
%
% Dependencies:
%     
R = [ -sin(lon0), -cos(lon0)*sin(lat0), cos(lon0)*cos(lat0); 
       cos(lon0), -sin(lon0)*sin(lat0), sin(lon0)*cos(lat0);
               0,            cos(lat0),           sin(lat0)];
        
%% TBD
X = [];  
Y = [];
Z = [];
end

