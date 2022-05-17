% haversin_distance.m
% 15/02/2013
function [dist dist_x dist_y] = haversin_distance( lat1, lat2, long1, long2, radius )
    dist = 2 * radius * asin( sqrt( ( sin( ( lat2 - lat1 ) / 2 ) ).^2 + cos( lat1 ) .* cos( lat2 ) .* ( sin( ( long2 - long1 ) / 2 ) ).^2 ) );
    dist_x = 2 * radius * asin( cos( lat1 ) .* sin( ( long2 - long1 ) / 2 ) );
    dist_y = 2 * radius * asin( sin( ( lat2 - lat1 ) / 2 ) );
end