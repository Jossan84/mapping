function [xEast,yNorth,zUp] = lla2enu(ellipsoid, lat, lon, h, lat0, lon0, h0)
% geodetic2enu
% 10/11/2021
% Transform geodetic coordinates to local east-north-up coordinates
%
% Usage:
%     [xEast,yNorth,zUp] = geodetic2enu(ellipsoid, lat, lon, h, lat0, lon0, h0)
%     arguments:
%         ellipsoid     Reference spheroid, only wgs84 and grs80 are implemented
%         lat      	    Geodetic latitude [rad]
%         lon           Geodetic longitude [rad]
%         h             Ellipsoidal height [m]
%         lat0      	  Reference geodetic latitude [rad]
%         lon0          Reference geodetic longitude [rad]
%         h0            Reference ellipsoidal height [m]
%
% Examples:
%     [xEast,yNorth,zUp] = lla2enu('wgs_84', lat, lon, h, lat0, lon0, h0)
%     [xEast,yNorth,zUp] = lla2enu('grs_80', lat, lon, h, lat0, lon0, h0)
%
% Dependencies:

R = nCe(lat0,lon0);
[xRef,yRef,zRef] = lla2ecef(ellipsoid, lat0, lon0, h0);

xEast = zeros(size(lat));
yNorth = zeros(size(lat));
zUp = zeros(size(lat));
for i = 1: length(lat)
    [xf,yf,zf] = lla2ecef(ellipsoid, lat(i), lon(i), h(i));
    X = R * [xf-xRef; yf-yRef; zf-zRef];

    xEast(i) = X(1); 
    yNorth(i) = X(2);
    zUp(i) = X(3);
end
end


function [Xe, Ye, Ze] = lla2ecef(ellipsoid, lat, lon, h)
% lla2ecef
% 11/11/2021
% Convert geodetic coordinates to Earth-centered Earth-fixed (ECEF) coordinates
%
% Usage:
%     [X,Y,Z] = geodetic2Ecef(ellipsoid, lat, lon, h)
%     arguments:
%         ellipsoid     Reference spheroid, only wgs84 and grs80 are implemented
%         lat      	    Geodetic latitude [rad]
%         lon           Geodetic longitude [rad]
%         h             Ellipsoidal height [m]
%
% Examples:
%     [X,Y,Z] = lla2ecef('wgs_84', lat, lon, h)
%     [X,Y,Z] = lla2ecef('grs_80', lat, lon, h)
%
% Dependencies:
%      
  
switch ellipsoid
  case 'wgs_84'
      a = 6378137;
      b = 6356752.31424518;
  case 'grs_80'
      a = 6378137.0;
      b = 6356752.314140; 
  otherwise
      error('Ellipsoid system is not implemented');
end
  
e = sqrt((a^2-b^2)/(a^2)); %[-] Excentricidad
    
N = (a/sqrt(1 - ((e^2) * sin(lat)^2)));
Xe = (h + N) * cos(lat) * cos(lon);
Ye = (h + N) * cos(lat) * sin(lon);    
Ze = ((((b^2)/(a^2)) * N) + h) * sin(lat);
end

% Rotation matrix form ECEF (Earth-Centered-Earth-Fixed) to 
% ENU (East-North-Up) reference frame.
function nCe = nCe(Lat,Lon)

    nCe = [         -sin(Lon),           cos(Lon),        0;
           -cos(Lon)*sin(Lat), -sin(Lon)*sin(Lat), cos(Lat);
            cos(Lon)*cos(Lat), sin(Lon)*cos(Lat),  sin(Lat)];
end