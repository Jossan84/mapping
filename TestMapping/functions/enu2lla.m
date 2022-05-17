function [lat, lon, h] = enu2lla(ellipsoid, xEast, yNorth, zUp, lat0, lon0, h0)
% enu2lla
% 11/11/2021
% Transform local east-north-up coordinates to geodetic coordinates
%
% Usage:
%     [lat, lon, h] = enu2lla(ellipsoid, xEast, yNorth, zUp, lat0, lon0, h0)
%     arguments:
%         ellipsoid    Reference spheroid, only wgs84 and grs80 are implemented
%         lat0      	 Reference geodetic latitude [rad]
%         lon0         Reference geodetic longitude [rad]
%         xEast        x ENU coordinates [m]
%         yNorth       y ENU coordinates [m]
%         zUp          z ENU coordinates [m]
%
% Examples:
%     [lat, lon, h] = enu2lla('wgs_84', xEast, yNorth, zUp, lat0, lon0, h0)
%
% Dependencies:
%     

error('enu2lla TBD');

               
function [X, Y, Z] = enu2ecef(ellipsoid, xEast, yNorth, zUp, lat0, lon0, h0)
% lla2ecef
% 11/11/2021
% Convert geodetic coordinates to Earth-centered Earth-fixed (ECEF) coordinates
%
% Usage:
%     [X,Y,Z] = geodetic2Ecef(ellipsoid, lat, lon, h)
%     arguments:
%         ellipsoid     Reference spheroid, only wgs84 and grs80 are implemented
%         lat0      	  Reference geodetic latitude [rad]
%         lon0          Reference geodetic longitude [rad]
%         h0            Reference ellipsoidal height [m]
%         xEast         x ENU coordinates [m]
%         yNorth        y ENU coordinates [m]
%         zUp           z ENU coordinates [m]
%
% Examples:
%     [X,Y,Z] = enu2ecef(wgs_84, xEast, yNorth, zUp, lat0, lon0, h0)
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

error('enu2ecef TBD');

end
               
% Rotation matrix form ENU (East-North-Up) reference frame to 
% ECEF (Earth-Centered-Earth-Fixed).
function nCe = nCe(Lat,Lon)
  nCe = [ -sin(lon0), -cos(lon0)*sin(lat0), cos(lon0)*cos(lat0); 
           cos(lon0), -sin(lon0)*sin(lat0), sin(lon0)*cos(lat0);
                   0,            cos(lat0),           sin(lat0)];
end               