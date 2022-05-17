function [ T ] = triangularWave( t, p, a)

    T = ((2*abs((t/p) - floor((t/p)+0.5))) * a);
    
end

