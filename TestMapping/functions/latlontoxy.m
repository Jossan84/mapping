function [x,y] = latlontoxy(reflat,reflon,lat,lon)
%latitude-longitude to xy conversion. 

re=6378137;
rp=6356800;
slat=sind(reflat*360/2);
clat=cosd(reflat);

r=(pi/180)*sqrt(re*re*clat*clat+rp*rp*slat*slat);
x=(lon-reflon)*r*clat;
y=(lat-reflat)*r;

 end