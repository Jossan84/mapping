function Y = firstOrderDynamicsOutputSum(y,u,n,a,b)
 
    Y = ((a^(n+1)-1)/(a-1) -1) * y + (1/(a-1))*((a^(n+1) - 1)/(a-1) - (n+1))* b * u;
    
    Y = Y - y*n;

end

