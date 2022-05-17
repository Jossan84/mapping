function [xEast,yNorth,zUp] = geodetic2Enu(ellipsoid, lat, lon, h, lat0, lon0, h0)
% geodetic2Enu
% 10/11/2021
% Transform geodetic coordinates to local east-north-up
%
% Usage:
%     [xEast,yEast,zEast] = geodetic2Enu(ellipsoid, lat, lon, h, lat0, lon0, h0)
%     arguments:
%         ellipsoid     Reference spheroid, only wgs84 and grs80 are implemented
%         lat      	    Geodetic latitude [rad]
%         lon           Geodetic longitude [rad]
%         h             Ellipsoidal height [m]
%         lat0      	Reference geodetic latitude [rad]
%         lon0          Reference geodetic longitude [rad]
%         h0            Reference ellipsoidal height [m]
%
% Examples:
%     [xEast,yEast,zEast] = geodetic2Enu('wgs_84', lat, lon, h, lat0, lon0, h0)
%     [xEast,yEast,zEast] = geodetic2Enu('grs_80', lat, lon, h, lat0, lon0, h0)
%
% Dependencies:
%     
% Step1: Convert geodetic coordinates to ECEF coordinates  
[X, Y, Z] = geodetic2Ecef(ellipsoid, lat, lon, h);

% Step 2: Convert ECEF coordinates to local ENU coordinates
[xEast,yNorth,zUp] = ecef2Enu(X, Y, Z, lat0, lon0);
  
end