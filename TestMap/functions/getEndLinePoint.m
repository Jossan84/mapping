function [endPointX, endPointY] = getEndLinePoint(startPointX, startPointY, endPointX, phi)

    m = tan(phi);
    n =  startPointY - (tan (phi) * startPointX );
    endPointY = m*endPointX + n;

end

