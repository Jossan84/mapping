function [s,t,n] = time2Reference(goalPoint,pointDistance,vx,T)

   x =  goalPoint(1);
   
   s = pointDistance * sqrt(1+((4/3)*(pointDistance - x)/(pointDistance + x)));
   
   t = s/vx;
   
   n = floor(t/T);

end

