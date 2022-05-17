function [clothoid1, clothoid2] = getClothoidCircleClothoidPoints(P1, P2, V, lambda1, lambda2, R, phi1_0, phi2_0)
    P1V = V - P1;
    P2V = V - P2;

%     phi1_0 = deg2rad(getAngle(P1, V));
%     phi2_0 = deg2rad(getAngle(P2, V));
    beta = acos( (P1V(1)*P2V(1) + P1V(2)*P2V(2)) / (sqrt(P1V(1)^2 + P1V(2)^2) * sqrt(P2V(1)^2 + P2V(2)^2)));
%     lambda1 = -1;
%     lambda2 = 1;
    omega = deg2rad(0);
    tau = (pi - (omega + beta))/2;
%     R = 1;
    L = 2*R*tau;
    A = sqrt(R*L);
    N = 100;
    s = linspace(0, L, N);
    ds = L/N;

    termXL = 0;
    termYL = 0;
    for k = 1 : N
        termXL = termXL + cos(s(k)^2 / (2*A^2));
        termYL = termYL + sin(s(k)^2 / (2*A^2));
    end

    XL = ds * termXL;
    YL = ds * termYL;

    distT1J1 = XL;
    distJ1H1 = YL * tan(tau);
    distJ2H2 = distJ1H1;
    distH1V = (R + YL/cos(tau)) * (sin(omega/2)/sin(beta/2));
    distH2V = distH1V;
    distT1V = distT1J1 + distJ1H1 + distH1V;
    distT2V = distT1V;

    T1 = V - distT1V * P1V/norm(P1V);
    T2 = V - distT2V * P2V/norm(P2V);

    clothoid1.x = zeros(length(s),1);
    clothoid1.y = zeros(length(s),1);
    clothoid1.x(1) = T1(1);
    clothoid1.y(1) = T1(2);
    clothoid2.x = zeros(length(s),1);
    clothoid2.y = zeros(length(s),1);
    clothoid2.x(1) = T2(1);
    clothoid2.y(1) = T2(2);


    for i = 1: length(s)
        clothoid1.x(i+1) = clothoid1.x(i) + ds * cos((lambda1*(s(i)^2)/(2*A^2)) + phi1_0);
        clothoid1.y(i+1) = clothoid1.y(i) + ds * sin((lambda1*(s(i)^2)/(2*A^2)) + phi1_0);

        clothoid2.x(i+1) = clothoid2.x(i) + ds * cos((lambda2*(s(i)^2)/(2*A^2)) + phi2_0);
        clothoid2.y(i+1) = clothoid2.y(i) + ds * sin((lambda2*(s(i)^2)/(2*A^2)) + phi2_0);
    end
    
end

