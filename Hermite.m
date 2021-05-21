function [xpositions, v, phi] = Hermite(route, n, k)
%Hermite används för att polynomanpassa y-värderna för x-värderna.
positions = 100*k;
xpositions = [];
v = [];
phi=[];
x1 = route.t(n);
x2 = route.t(n + 1);


for i = 1:positions
    x = x1 + (x2-x1)*i/positions;
    p = herm(    x,...
                    [x1, x2],...
                    [route.r(n), route.r(n + 1)],...
                    [route.rd(n), route.rd(n + 1)]);
    p2= herm(    x,...
                    [x1, x2],...
                    [route.phi(n), route.phi(n + 1)],...
                    [route.phid(n), route.phid(n + 1)]);         

    v = [v p];
    phi=[phi,p2];
    xpositions = [xpositions x];
end

function yvec = herm(xq,x,y,k)

x1=x(1);
x2=x(2);
y1=y(1);
y2=y(2);
k1=k(1);
k2=k(2);

h = x2-x1;

c1 = y1 ;
c2 = (y2-y1)/h;
c3 = (k2-c2)/h^2;
c4 = (k1-c2)/h^2;

yvec = c1+c2*(xq-x1)+c3*(xq-x1)^2*(xq-x2)+c4*(xq-x1)*(xq-x2)^2;
        
end
end