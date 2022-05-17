function [path] = Path(varargin)
            
    for i = 1 :2: nargin
        switch varargin{i}
            case 'Bezier'
                degree = '1';
                n = 3;
                p0 = [ 0 0   100; 0   100 100];   % [x;y] [m];
                p1 = [ 0 100 100; 100 100 200];   %  " "   "
%                 p0 = [ 0 ; 0  ];   % [x;y] [m];
%                 p1 = [ 0 ; 100];   %  " "   "
                p2 = [];
                p3 = [];
                p4 = [];
                p5 = [];
            case 'Bezier2'
                degree = '2';
                n = 3;
                p0 = [ 0 0   100; 0   100 100];   % [x;y] [m];
                p1 = [ 0 50  100; 50  100 150];   %  " "   " 
                p2 = [ 0 100 100; 100 100 200];   %  " "   "
%                 p0 = [ 0 ; 0  ];   % [x;y] [m];
%                 p1 = [ 0 ; 50 ];   %  " "   " 
%                 p2 = [ 0 ; 100];   %  " "   "
                p3 = [];
                p4 = [];
                p5 = [];
            case 'Bezier3'
                degree = '3';
                n = 3;
                p0 = [ 0 0   100; 0   100 100];   % [x;y] [m];
                p1 = [ 0 25  100; 25  100 125];   %  " "   " 
                p2 = [ 0 75  100; 75  100 175];   %  " "   " 
                p3 = [ 0 100 100; 100 100 200];   %  " "   "
%                 p0 = [ 0 ; 0  ];   % [x;y] [m];
%                 p1 = [ 0 ; 25 ];   %  " "   " 
%                 p2 = [ 0 ; 75 ];   %  " "   " 
%                 p3 = [ 0 ; 100];   %  " "   "
                p4 = [];
                p5 = [];
            case 'Bezier5'
                degree = '5';
                n = 3;
                p0 = [ 0 0   100; 0   100 100];   % [x;y] [m];
                p1 = [ 0 20  100; 20  100 120];   %  " "   " 
                p2 = [ 0 40  100; 40  100 140];   %  " "   " 
                p3 = [ 0 60  100; 60  100 160];   %  " "   " 
                p4 = [ 0 80  100; 80  100 180];   %  " "   " 
                p5 = [ 0 100 100; 100 100 200];   %  " "   "
%                 p0 = [ 0 ; 0  ];   % [x;y] [m];
%                 p1 = [ 0 ; 20 ];   %  " "   " 
%                 p2 = [ 0 ; 40 ];   %  " "   " 
%                 p3 = [ 0 ; 60 ];   %  " "   " 
%                 p4 = [ 0 ; 80 ];   %  " "   " 
%                 p5 = [ 0 ; 100];   %  " "   "
            otherwise
                error('Wrong argument');
        end
    end
    
    points = [p0 p1 p2 p3 p4 p5];
    
    path.points                 = points;
    path.degree                 = degree;
    path.n                      = n;
    
    path.getPath                = @getPath;
    path.getCrossTrackError     = @getCrossTrackError;
    path.getLookAheadPoint      = @getLookAheadPoint;
    path.plotPath               = @plotPath;
    
end

function [x, y] = getPath(path,t,segment)

        switch path.degree
            case '1'
                P0 = path.points(:,segment);
                P1 = path.points(:,segment+path.n);
                x  = @(t) (1-t)*P0(1) + t*P1(1);
                y  = @(t) (1-t)*P0(2) + t*P1(2);
                x  = x(t);
                y  = y(t);                
                
            case '2'
                P0 = path.points(:,segment);
                P1 = path.points(:,segment+path.n);
                P2 = path.points(:,segment+path.n*2);
                x = @(t) ((1-t).^2)*P0(1) + 2*(1-t).*t*P1(1) + (t.^2)*P2(1);
                y = @(t) ((1-t).^2)*P0(2) + 2*(1-t).*t*P1(2) + (t.^2)*P2(2);
                x = x(t);
                y = y(t);
            case '3'
                P0 = path.points(:,segment);
                P1 = path.points(:,segment+path.n);
                P2 = path.points(:,segment+path.n*2);
                P3 = path.points(:,segment+path.n*3);
                x = @(t) ((1-t).^3)*P0(1) + 3*((1-t).^2).*t*P1(1) + 3*(1-t).*...
                         (t.^2)*P2(1) + (t.^3)*P3(1);
                y = @(t) ((1-t).^3)*P0(2) + 3*((1-t).^2).*t*P1(2) + 3*(1-t).*...
                         (t.^2)*P2(2) + (t.^3)*P3(2);
                x = x(t);
                y = y(t);
            case '5'
                P0 = path.points(:,segment);
                P1 = path.points(:,segment+path.n);
                P2 = path.points(:,segment+path.n*2);
                P3 = path.points(:,segment+path.n*3);
                P4 = path.points(:,segment+path.n*4);
                P5 = path.points(:,segment+path.n*5);
                x = @(t) ((1-t).^5)*P0(1) + 5*((1-t).^4).*t*P1(1) + 10*((1-t).^3).*...
                         (t.^2)*P2(1) + 10*((1-t).^2).*(t.^3)*P3(1) + 5*(1-t).*...
                         (t.^4)*P4(1) + (t.^5)*P5(1);
                y = @(t) ((1-t).^5)*P0(2) + 5*((1-t).^4).*t*P1(2) + 10*((1-t).^3).*...
                         (t.^2)*P2(2) + 10*((1-t).^2).*(t.^3)*P3(2) + 5*(1-t).*...
                         (t.^4)*P4(2) + (t.^5)*P5(2);
                x = x(t);
                y = y(t);
            otherwise
                error('Wrong argument');
        end
end

function [x,y] = plotPath(path,nSteps)

    t      = linspace(0,1,nSteps);
    
    x = [];
    y = [];
    
    for i=1:path.n
        [x_i, y_i] = path.getPath(path,t,i);       
        x = [x x_i];
        y = [y y_i];
    end    
    
    figure(1);
    plot(x,y,'b.-');
    title('Path');
    ylabel('y [m]');
    xlabel('x [m]');
    
end
function [cte,segment, e, nIter] = getCrossTrackError(path,x0,y0,method)

        switch path.degree
            case '1'
                for i=1:path.n
                    P0 = path.points(:,i);
                    P1 = path.points(:,i+path.n);
                
                    x   = @(t) ((1-t)*P0(1) + t*P1(1));
                    dx  = @(t) (P1(1) - P0(1));
                    ddx = @(t) (0);
                
                    y   = @(t) ((1-t)*P0(2) + t*P1(2));
                    dy  = @(t) (P1(2) - P0(2));
                    ddy = @(t) (0);
                    
                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2));
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t) +...
                               (y(t)-y0)*ddy(t)) );
        
                    t0 = 0.5;
                    errMax = 0.01;
        
                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f1,f2);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end
                        cte_i(i) = sqrt(f(t));  
                    
                end
                        [cte, segment] = min(cte_i);
            case '2'
                for i=1:path.n
                    P0  = path.points(:,i);
                    P1  = path.points(:,i+path.n);
                    P2  = path.points(:,i+path.n*2);

                    x   = @(t) ((1-t).^2)*P0(1) + 2*(1-t).*t*P1(1) + (t.^2)*P2(1);
                    dx  = @(t) (2*P2(1).*t - 2*P1(1).*t + P0(1)*(2.*t - 2) - P1(1)*...
                               (2.*t - 2));
                    ddx = @(t) (2*P0(1) - 4*P1(1) + 2*P2(1));

                    y   = @(t) ((1-t).^2)*P0(2) + 2*(1-t).*t*P1(2) + (t.^2)*P2(2);
                    dy  = @(t) (2*P2(2).*t - 2*P1(2).*t + P0(2)*(2.*t - 2) - P1(2)*...
                               (2.*t - 2));
                    ddy = @(t) (2*P0(2) - 4*P1(2) + 2*P2(2));

                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2));
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t) +...
                               (y(t)-y0)*ddy(t)) );
        
                    t0 = 0.5;
                    errMax = 0.01;
        
                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f1,f2);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end
                        cte_i(i) = sqrt(f(t)); 
                end
                        [cte, segment] = min(cte_i);
            case '3'
                for i=1:path.n
                    P0  = path.points(:,i);
                    P1  = path.points(:,i+path.n);
                    P2  = path.points(:,i+path.n*2);
                    P3  = path.points(:,i+path.n*3);

                    x   = @(t) ((1-t).^3)*P0(1) + 3*((1-t).^2).*t*P1(1) + 3*(1-t).*...
                               (t.^2)*P2(1)...
                                + (t.^3)*P3(1);
                    dx  = @(t) (3*P3(1).*t.^2 - 3*P2(1).*t.^2 - 3*P0(1).*(t - 1).^2 ...
                               + 3*P1(1).*...
                               (t - 1).^2 + 3*P1(1).*t*(2*t - 2) - 2*P2(1)*t*(3*t - 3));
                    ddx = @(t) (6*P1(1).*t - 12*P2(1).*t + 6*P3(1).*t - 3*P0(1)*...
                               (2.*t - 2) + ...
                               6*P1(1)*(2.*t - 2) - 2*P2(1)*(3*t - 3));

                    y   = @(t) ((1-t).^3)*P0(2) + 3*((1-t).^2).*t*P1(2) + 3*(1-t).*...
                               (t.^2)*P2(2) +...
                               (t.^3)*P3(2);
                    dy  = @(t) (3*P3(2).*t.^2 - 3*P2(2).*t.^2 - 3*P0(2).*(t - 1).^2 ...
                               + 3*P1(2).*...
                               (t - 1).^2 + 3*P1(2).*t*(2*t - 2) - 2*P2(2)*t*(3*t - 3));
                    ddy = @(t) (6*P1(2).*t - 12*P2(2).*t + 6*P3(2).*t - 3*P0(2)*...
                               (2.*t - 2) +...
                               6*P1(2)*(2.*t - 2) - 2*P2(2)*(3*t - 3));
                           
                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2));
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t) +...
                               (y(t)-y0)*ddy(t)) );
        
                    t0 = 0.5;
                    errMax = 0.01;
        
                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f1,f2);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end
                        cte_i(i) = sqrt(f(t)); 
                end
                        [cte, segment] = min(cte_i);
                
            case '5'
                for i=1:path.n
                    P0  = path.points(:,i);
                    P1  = path.points(:,i+path.n);
                    P2  = path.points(:,i+path.n*2);
                    P3  = path.points(:,i+path.n*3);
                    P4  = path.points(:,i+path.n*4);
                    P5  = path.points(:,i+path.n*5);
                    x   = @(t) ((1-t).^5)*P0(1) + 5*((1-t).^4).*t*P1(1) + 10*((1-t).^3).*...
                               (t.^2)...
                               *P2(1) + 10*((1-t).^2).*(t.^3)*P3(1) + 5*(1-t).*(t.^4)*P4(1)+...
                               (t.^5)*P5(1);
                    dx  = @(t) (5*P5(1)*t^4 - 5*P4(1)*t^4 - 5*P0(1)*(t - 1)^4 + 5*P1(1)*(t - 1)^4 ...
                               + 20*P1(1)*t*(t - 1)^3 - 20*P2(1)*t*(t - 1)^3 + 10*P3(1)*t^3* ...
                               (2*t - 2) - 4*P4(1)*t^3*(5*t - 5)- 30*P2(1)*t^2*(t - 1)^2 + ...
                               30*P3(1)*t^2*(t - 1)^2);
                    ddx = @(t) (5*P5(1)*t^4 - 5*P4(1)*t^4 - 5*P0(1)*(t - 1)^4 + 5*P1(1)*(t - 1)^4 ...
                               + 20*P1(1)*t*(t - 1)^3 - 20*P2(1)*t*(t - 1)^3 + 10*P3(1)*t^3* ...
                               (2*t - 2) - 4*P4(1)*t^3*(5*t - 5) - 30*P2(1)*t^2*(t - 1)^2 + ...
                               30*P3(1)*t^2*(t - 1)^2);       


                    y   = @(t) ((1-t).^5)*P0(2) + 5*((1-t).^4).*t*P1(2) + 10*((1-t).^3).*(t.^2)*P2(2)...
                               + 10*((1-t).^2).*(t.^3)*P3(2) + 5*(1-t).*(t.^4)*P4(2) + (t.^5)*P5(2);
                    dy  = @(t) (5*P5(2)*t^4 - 5*P4(2)*t^4 - 5*P0(2)*(t - 1)^4 + 5*P1(2)*(t - 1)^4 ...
                               + 20*P1(2)*t*(t - 1)^3 - 20*P2(2)*t*(t - 1)^3 + 10*P3(2)*t^3* ...
                               (2*t - 2) - 4*P4(2)*t^3*(5*t - 5)- 30*P2(2)*t^2*(t - 1)^2 + ...
                               30*P3(2)*t^2*(t - 1)^2);
                    ddy = @(t) (5*P5(2)*t^4 - 5*P4(2)*t^4 - 5*P0(2)*(t - 1)^4 + 5*P1(2)*(t - 1)^4 ...
                               + 20*P1(2)*t*(t - 1)^3 - 20*P2(2)*t*(t - 1)^3 + 10*P3(2)*t^3* ...
                               (2*t - 2) - 4*P4(2)*t^3*(5*t - 5) - 30*P2(2)*t^2*(t - 1)^2 ...
                               + 30*P3(2)*t^2*(t - 1)^2);
                           
                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2));
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t) +...
                               (y(t)-y0)*ddy(t)) );
        
                    t0 = 0.5;
                    errMax = 0.01;
        
                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f1,f2);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end
                        cte_i(i) = sqrt(f(t));
                end
                        [cte, segment] = min(cte_i);
            otherwise
                error('Wrong argument');
        end
                       
end

function [goalPoint, e, nIter] = getLookAheadPoint(path,x0,y0,R,method)

        switch path.degree
            case '1'
                for i=1:1
                    P0 = path.points(:,i);
                    P1 = path.points(:,i+path.n);
                
                    x   = @(t) ((1-t)*P0(1) + t*P1(1));
                    dx  = @(t) (P1(1) - P0(1));
                    ddx = @(t) (0);

                    y   = @(t) ((1-t)*P0(2) + t*P1(2));
                    dy  = @(t) (P1(2) - P0(2));
                    ddy = @(t) (0);
                    
                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2) - R^2);
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t)...
                               + (y(t)-y0)*ddy(t)) );
               
                    t0 = 0.5;
                    errMax = 0.01;  % [%]

                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f,f1);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end

                    goalPoint(1) = x(t);
                    goalPoint(2) = y(t);
                end
                
            case '2'
                for i=1:1
                    P0  = path.points(:,i);
                    P1  = path.points(:,i+path.n);
                    P2  = path.points(:,i+path.n*2);

                    x   = @(t) ((1-t).^2)*P0(1) + 2*(1-t).*t*P1(1) + (t.^2)*P2(1);
                    dx  = @(t) (2*P2(1).*t - 2*P1(1).*t + P0(1)*(2.*t - 2) - P1(1)*...
                               (2.*t - 2));
                    ddx = @(t) (2*P0(1) - 4*P1(1) + 2*P2(1));

                    y   = @(t) ((1-t).^2)*P0(2) + 2*(1-t).*t*P1(2) + (t.^2)*P2(2);
                    dy  = @(t) (2*P2(2).*t - 2*P1(2).*t + P0(2)*(2.*t - 2) - P1(2)*...
                               (2.*t - 2));
                    ddy = @(t) (2*P0(2) - 4*P1(2) + 2*P2(2));
                    
                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2) - R^2);
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t)...
                               + (y(t)-y0)*ddy(t)) );
               
                    t0 = 0.5;
                    errMax = 0.01;  % [%]

                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f,f1);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end

                    goalPoint(1) = x(t);
                    goalPoint(2) = y(t);                    
                    
                end 
            case '3'
                for i=1:1
                    P0  = path.points(:,i);
                    P1  = path.points(:,i+path.n);
                    P2  = path.points(:,i+path.n*2);
                    P3  = path.points(:,i+path.n*3);

                    x   = @(t) ((1-t).^3)*P0(1) + 3*((1-t).^2).*t*P1(1) + 3*(1-t).*...
                               (t.^2)*P2(1)...
                                + (t.^3)*P3(1);
                    dx  = @(t) (3*P3(1).*t.^2 - 3*P2(1).*t.^2 - 3*P0(1).*(t - 1).^2 ...
                               + 3*P1(1).*...
                               (t - 1).^2 + 3*P1(1).*t*(2*t - 2) - 2*P2(1)*t*(3*t - 3));
                    ddx = @(t) (6*P1(1).*t - 12*P2(1).*t + 6*P3(1).*t - 3*P0(1)*...
                               (2.*t - 2) + ...
                               6*P1(1)*(2.*t - 2) - 2*P2(1)*(3*t - 3));

                    y   = @(t) ((1-t).^3)*P0(2) + 3*((1-t).^2).*t*P1(2) + 3*(1-t).*...
                               (t.^2)*P2(2) +...
                               (t.^3)*P3(2);
                    dy  = @(t) (3*P3(2).*t.^2 - 3*P2(2).*t.^2 - 3*P0(2).*(t - 1).^2 ...
                               + 3*P1(2).*...
                               (t - 1).^2 + 3*P1(2).*t*(2*t - 2) - 2*P2(2)*t*(3*t - 3));
                    ddy = @(t) (6*P1(2).*t - 12*P2(2).*t + 6*P3(2).*t - 3*P0(2)*...
                               (2.*t - 2) +...
                               6*P1(2)*(2.*t - 2) - 2*P2(2)*(3*t - 3));
                           
                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2) - R^2);
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t)...
                               + (y(t)-y0)*ddy(t)) );
               
                    t0 = 0.5;
                    errMax = 0.01;  % [%]

                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f,f1);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end

                    goalPoint(1) = x(t);
                    goalPoint(2) = y(t);                            
                           
                end
            case '5'
                for i=1:1
                    P0  = path.points(:,i);
                    P1  = path.points(:,i+path.n);
                    P2  = path.points(:,i+path.n*2);
                    P3  = path.points(:,i+path.n*3);
                    P4  = path.points(:,i+path.n*4);
                    P5  = path.points(:,i+path.n*5);
                    
                    x   = @(t) ((1-t).^5)*P0(1) + 5*((1-t).^4).*t*P1(1) + 10*((1-t).^3).*...
                               (t.^2)*P2(1) + 10*((1-t).^2).*(t.^3)*P3(1) + 5*(1-t).*...
                               (t.^4)*P4(1) + (t.^5)*P5(1);
                    dx  = @(t) (5*P5(1)*t^4 - 5*P4(1)*t^4 - 5*P0(1)*(t - 1)^4 + 5*P1(1)*(t - 1)^4 ...
                               + 20*P1(1)*t*(t - 1)^3 - 20*P2(1)*t*(t - 1)^3 + 10*P3(1)*t^3* ...
                               (2*t - 2) - 4*P4(1)*t^3*(5*t - 5)- 30*P2(1)*t^2*(t - 1)^2 + ...
                               30*P3(1)*t^2*(t - 1)^2);
                    ddx = @(t) (5*P5(1)*t^4 - 5*P4(1)*t^4 - 5*P0(1)*(t - 1)^4 + 5*P1(1)*(t - 1)^4 ...
                               + 20*P1(1)*t*(t - 1)^3 - 20*P2(1)*t*(t - 1)^3 + 10*P3(1)*t^3* ...
                               (2*t - 2) - 4*P4(1)*t^3*(5*t - 5) - 30*P2(1)*t^2*(t - 1)^2 + ...
                               30*P3(1)*t^2*(t - 1)^2);       


                    y   = @(t) ((1-t).^5)*P0(2) + 5*((1-t).^4).*t*P1(2) + 10*((1-t).^3).*(t.^2)*P2(2)...
                               + 10*((1-t).^2).*(t.^3)*P3(2) + 5*(1-t).*(t.^4)*P4(2) + (t.^5)*P5(2);
                    dy  = @(t) (5*P5(2)*t^4 - 5*P4(2)*t^4 - 5*P0(2)*(t - 1)^4 + 5*P1(2)*(t - 1)^4 ...
                               + 20*P1(2)*t*(t - 1)^3 - 20*P2(2)*t*(t - 1)^3 + 10*P3(2)*t^3* ...
                               (2*t - 2) - 4*P4(2)*t^3*(5*t - 5)- 30*P2(2)*t^2*(t - 1)^2 + ...
                               30*P3(2)*t^2*(t - 1)^2);
                    ddy = @(t) (5*P5(2)*t^4 - 5*P4(2)*t^4 - 5*P0(2)*(t - 1)^4 + 5*P1(2)*(t - 1)^4 ...
                               + 20*P1(2)*t*(t - 1)^3 - 20*P2(2)*t*(t - 1)^3 + 10*P3(2)*t^3* ...
                               (2*t - 2) - 4*P4(2)*t^3*(5*t - 5) - 30*P2(2)*t^2*(t - 1)^2 ...
                               + 30*P3(2)*t^2*(t - 1)^2);

                    % Objective function and derivatives to minimization
                    f   = @(t) ( ((x(t)-x0)^2) + ((y(t)-y0)^2) - R^2);
                    f1  = @(t) ( 2*(x(t)-x0) * dx(t)  + 2*(y(t)-y0) * dy(t) );
                    f2  = @(t) ( 2*(dx(t)*dx(t) + (x(t)-x0)*ddx(t) ) + 2*(dy(t)*dy(t)...
                               + (y(t)-y0)*ddy(t)) );
               
                    t0 = 0.5;
                    errMax = 0.01;  % [%]

                    % Compute minimization
                    switch method
                        case 'Newton'
                            [t, e, nIter] = newtonMethod(t0,errMax ,f,f,f1);
                        case 'Bisection'
                            [t, e, nIter] = bisectionMethod(errMax, f);
                        case 'RegulaFalsi'
                            [t, e, nIter] = regulaFalsiMethod(errMax, f);
                        otherwise
                            error('Method not valid');   
                    end

                    goalPoint(1) = x(t);
                    goalPoint(2) = y(t); 
                end
            otherwise
                error('Wrong argument');
        end
                       
end
