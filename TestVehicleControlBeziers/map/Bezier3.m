function [bezier3] = Bezier3(varargin)
    
    % Default params
    params.nPoints       = 100;       % []
    params.p1            = [ 5;-10];  % [x;y] [m];
    params.p2            = [10; 18];  %  " "   "
    params.p3            = [30; 40];  %  " "   "
    params.p4            = [45; 15];  %  " "   " 
        
    for i = 1 :2: nargin
        switch varargin{i}
            case 'nPoints'
                params.nPoints = varargin{i + 1};
            case 'p1'
                params.p1 = varargin{i + 1}; 
            case 'p2'
                params.p2 = varargin{i + 1};
            case 'p3'
                params.p3 = varargin{i + 1};
            case 'p4'
                params.p4 = varargin{i + 1};
            otherwise
                error('Wrong argument');
        end
    end
    
    params.t  = linspace(0,1,params.nPoints);
    
    [~,m] = size(params.t);
    
    points    = zeros(2,m);
    deltaYaw  = zeros(1,m);
    curvature = zeros(1,m);

    bezier3.params     = params;
    bezier3.points     = points;
    bezier3.deltaYaw   = deltaYaw;
    bezier3.curvature  = curvature;
    
    bezier3.getPoints        = @getPoints;
    bezier3.getDeltaYaw      = @getDeltaYaw;
    bezier3.getCurvature     = @getCurvature;
    bezier3.plotBezier       = @plotBezier;
    
end

function [bezier3] = getPoints(bezier3)

    a = (1-bezier3.params.t).^3;
    b = 3*(1-bezier3.params.t).^2.*bezier3.params.t;
    c = 3*(1-bezier3.params.t).*bezier3.params.t.^2;
    d = bezier3.params.t.^3;
 
    [~,m] = size(a);

    P1 = bezier3.params.p1 * ones([1 m]); 
    P2 = bezier3.params.p2 * ones([1 m]);
    P3 = bezier3.params.p3 * ones([1 m]);
    P4 = bezier3.params.p4 * ones([1 m]);
 
    A  = [a; a];
    B  = [b; b];
    C  = [c; c];
    D  = [d; d];
    
    % Get Bezier points whith Kronecker Product 
    bezier3.points = (A.*P1) + (B.*P2) + (C.*P3) + (D.*P4);   

end

function [bezier3] = getDeltaYaw(bezier3)

% ************************* TBD ************************************
     disp('TBD');
%     [~,m] = size(bezier3.params.t);
%     
%     a = (-3*bezier3.params.t.^2) + 6*bezier3.params.t -3*ones(1,m);
%     b = (9*bezier3.params.t.^2) - 12*bezier3.params.t +3*ones(1,m);
%     c = (-9*bezier3.params.t.^2) + 6*bezier3.params.t;
%     d = (3*bezier3.params.t.^2);
%  
%     P1 = bezier3.params.p1 * ones([1 m]);
%     P2 = bezier3.params.p2 * ones([1 m]);
%     P3 = bezier3.params.p3 * ones([1 m]);
%     P4 = bezier3.params.p4 * ones([1 m]);
%  
%     A  = [a; a];
%     B  = [b; b];
%     C  = [c; c];
%     D  = [d; d];
%     
%     % Get Bezier points whith Kronecker Product 
%     bezier3.deltaYaw = (A.*P1) + (B.*P2) + (C.*P3) + (D.*P4);
    
end

function [bezier3] = getCurvature(bezier3)
% ************************* TBD ************************************
    disp('TBD');
end

function [] = plotBezier(bezier3)

%% Plot Bezier
hold on

  plot(bezier3.points(1,:),bezier3.points(2,:),'b.-');
  text(bezier3.params.p1(1),bezier3.params.p1(2),'p_1');
  text(bezier3.params.p2(1),bezier3.params.p2(2),'p_2');
  text(bezier3.params.p3(1),bezier3.params.p3(2),'p_3');
  text(bezier3.params.p4(1),bezier3.params.p4(2),'p_4');
  
  plot(bezier3.params.p1(1),bezier3.params.p1(2),'k.');
  plot(bezier3.params.p2(1),bezier3.params.p2(2),'k.');
  plot(bezier3.params.p3(1),bezier3.params.p3(2),'k.');
  plot(bezier3.params.p4(1),bezier3.params.p4(2),'k.');
  
  %xlim([0 50]);
  %ylim([-20 50]);
  axis equal;

hold off
    
end
