% Funktion för RK4
function [route] = RK4 (h, H)
u = [H+1 0 0 0];                          % B.V
t = 0;                                    % Starttid
route = struct(    't',        t,...  %Structure array
                        'r',        u(1),...
                        'rd',     u(2),...
                        'phi',      u(3),...
                        'phid',   u(4));
%Kollar när r' är negativt: Raketen minskar i höjd.
%RK4
while route.rd <= 0                        
    k1 = h*Function(t, u, H);
    k2 = h*Function(t+0.5*h, u+0.5*k1, H);
    k3 = h*Function(t+0.5*h, u+0.5*k2, H);
    k4 = h*Function(t+h, u+k3, H);

    u = u + (k1+2*k2+2*k3+k4)/6;
    t = t + h;   
    route.t = [route.t; t];
    route.r = [route.r; u(1)];
    route.rd = [route.rd; u(2)];
    route.phi = [route.phi; u(3)];
    route.phid = [route.phid; u(4)];  
end