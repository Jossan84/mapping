function [point] = getIntersectionPoint(line1, line2)

    x = ((line1.x1*line1.y2 - line1.y1*line1.x2)*(line2.x1-line2.x2) - (line1.x1-line1.x2) * (line2.x1*line2.y2 - line2.y1*line2.x2)) / ...
        ((line1.x1-line1.x2) * (line2.y1-line2.y2) - (line1.y1-line1.y2) * (line2.x1-line2.x2));
    y = ((line1.x1*line1.y2 - line1.y1*line1.x2)*(line2.y1-line2.y2) - (line1.y1-line1.y2) * (line2.x1*line2.y2 - line2.y1*line2.x2)) / ...
        ((line1.x1-line1.x2) * (line2.y1-line2.y2) - (line1.y1-line1.y2) * (line2.x1-line2.x2));
    point = [x, y];
    
end

