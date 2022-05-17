function [ t, e, nIter] = bisectionMethod( err, f )

    a = 0;
    b = 1;
    t = 0;
    nIter = 0;
    e = 100;
    while e>err
    
        ti = 0.5*(a + b);
        
        if ((f(a)*f(ti)) < 0)
            a = a;
            b = ti;
        else
            a = ti;
            b = b;        
        end
        
       e     = abs(((ti-t)/t)*100); % Error [%]
       t     = ti;
        
    nIter = nIter + 1;              % Nº iterations            
    end
    

end

