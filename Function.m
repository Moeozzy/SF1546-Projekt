%Uppgift a)

function [u_vec] = Function(t, u, H)


global alpha 
%Angivna värden
g = 20;                         
F = (g)/(1+H)^2;  
%Ekvation skriven till ett system av första ordningens
%begynnelsevärdesproblem.
x1 = u(2);
x2 = F*cosd(alpha)-g/(u(1)^2) + u(1)*(u(4))^2;
x3 = u(4);
x4 = (F*sind(alpha)-2*u(2)*u(4))/u(1);

u_vec = [x1 x2 x3 x4];