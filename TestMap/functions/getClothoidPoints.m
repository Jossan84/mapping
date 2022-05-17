function [x, y] = getClothoidPoints(F, C, lambda, phi0, alpha)

    R = sqrt((F(1)-C(1))^2 + ((F(2)-C(2))^2));
    
    % Clothoid
    if (lambda * (alpha-phi0)) >= 0
       tau = lambda * (alpha - phi0); 
    else
       tau = 2*pi + lambda * (alpha -phi0); 
    end
    L = 2 * R * tau;
    A = sqrt(R*L);
    N = 100;
    s = linspace(0, L, N);
    ds = ones(1,N) * L/N;

    % Note: If s not increase linearly, we can supose that we not travel at the
    %       same speed during all the clothoid. "TRY IN A TEST"
    % TEST
    % s1 = linspace(0,  L/3, N/2);
    % s2 = linspace(L/3,  L, N/2);
    % s = [s1, s2(2:end)];
    % ds = [diff(s)];
    % figure(1);
    % plot(s, 'b.-');
    % END TEST

    x = zeros(length(s),1);
    y = zeros(length(s),1);
    x(end) = F(1);
    y(end) = F(2);
    for i = length(s)-1 : -1: 1
        x(i) = x(i+1) - ds(i) * cos((lambda * (s(i)^2)/(2*A^2)) + phi0);
        y(i) = y(i+1) - ds(i) * sin((lambda * (s(i)^2)/(2*A^2)) + phi0);
    end

end

