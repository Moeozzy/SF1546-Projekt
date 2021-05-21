function [route] = Paths(route)
% Denna funktion beräknar de olika banorna för raketen genom parametrarna r,phi,t och v. genom Linearpol funktionen och Hermite för att answeepa. Detta görs när r' är lika med noll. 

% Feluppskattning för metoderna:
RK_r_err = 1.677e-06;
RK_phi_err = 7.7798e-07;
RK_phid_err = 2.0611e-06;
herm_rel_err = 0.034056e-6;

rd_sweep=0;

[t_sweep, t_err] = Linearpol(       route.rd(end-2:end),...
                                route.t(end-2:end),...
                                rd_sweep);
r_sweep=herm(    t_sweep,...
                route.t(end-1:end),...
                route.r(end-1:end),...
                route.rd(end-1:end));
phi_sweep=herm(  t_sweep,...
                route.t(end-1:end),...
                route.phi(end-1:end),...
                route.phid(end-1:end));
[phid_sweep, phid_err] =Linearpol(  route.rd(end-2:end),...
                                    route.phid(end-2:end),...
                                    rd_sweep);
        
v_sweep=phid_sweep*r_sweep;
v_sweep_err = phid_err + phid_sweep*RK_phid_err;
        
%Spara för att sedan föra över till tabeller.
route.t_sweep = t_sweep;
route.t_err = t_err;
route.r_sweep = r_sweep;
route.r_err = r_sweep*(herm_rel_err + RK_r_err);
route.phi_sweep = phi_sweep;
route.phi_err = phi_sweep*(herm_rel_err + RK_phi_err);
route.v_sweep = v_sweep;
route.v_err = v_sweep_err;
function yvec = herm(xq,x,y,k)
% Här används herm för att anpassa.
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