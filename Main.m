% Projektarbete "Raketen"
% Axel Arvidzon och Marcus Bohl
close all; clear all; clc;

global alpha

plus_minus = char(177);   % ±
h = 0.01;  % Steglängd
bisec_err = 0.0001;

% Upggfit b)

alpha = 90;
t_list = [];                        
r_list = [];                        
phi_list = [];                      
H_list = [];                        
crash_list = [];                    
startH = [5, 4, 3, 2];

% Plotta banor för olika H

figure (1)
plotStyle = {'b','g','r','m'};                       
for i = 1:length(startH)                   
    H = startH(i);
    trajectory = RK4(h, H);                       
    trajectory.H=H;
    routes(i) = trajectory;
    %Plotta banorna
    polar(trajectory.phi, trajectory.r, plotStyle{i})                      
    view([90 -90])                    
    grid on; hold all;
    legendInfo{i} = ['Start höjd ' num2str(H) ' Jordradie'];       
end
title ('banor för olika starthöjder, \alpha =90')
legend(legendInfo)     

%Hitta alla värden för banorna
for i = 1:length(startH)
    routes_temp(i)=Paths(routes(i));
end
routes=routes_temp;


% Svar på uppgift c, tabell med alla värden på de olika banornas lägsta
% punkt.
Table(routes)

% Uppgift d) Hitta H* för när raketen precis passerar jorden.
% Bisektionsmetoden används för att hitta H*

start1=1;       %Start gissning
start2=4;        %Slut gissning
start=[start1 start2];

H_critical = Bisection(@Diff_Alpha, start, h);
% Beräkna banan för H*
route_critical = RK4(h, H_critical);
route_critical = Paths(route_critical);
fprintf('Raketen passerar precis jordytan vid: %0.3f%c%f jordradier\n', H_critical, plus_minus, bisec_err)
fprintf('passeringsradie: %0.3f%c%f jordradier\n', route_critical.r_sweep, plus_minus, route_critical.r_err)
fprintf('hastigheten: %0.3f%c%f jordradier/h\n', route_critical.v_sweep, plus_minus, route_critical.v_err)
% Banlängden
Y.critical=route_critical.r.*sin(route_critical.phi);
X.critical=route_critical.r.*cos(route_critical.phi); %Polära koordinater till kartetiska
traj_length = Path_Lenght(X.critical, 0, Y.critical); %Avstånd
traj_length_2h = Path_Lenght((X.critical(1:2:end)), 0, Y.critical(1:2:end));
% Svaret på uppgift e)
fprintf('Banlängd: %0.3f%c%f jordradier\n', traj_length, plus_minus, abs(traj_length - traj_length_2h))

%Plotten för H*
figure(2)
polar(route_critical.phi, route_critical.r,'r')
title('Bana för H*, \alpha =90')
view([90 -90])
hold on;
%Plot för jorden
phi_earth = 0:360/(length(route_critical.r)):360;
r_earth = ones(1,length(route_critical.r)+1);
polar(phi_earth,r_earth,'b')
legend('Bana','Jorden')
v0 = route_critical.v_sweep;
% Plotten för bankuurvan
figure(3) 
plot (X.critical, Y.critical)
leg = {'RK4'};
legend(leg)
title('Interpolerad kurva för raketen mot kurva från RK4')
xlabel('x [jordradier]')
ylabel('y [jordradier]')

% Uppgift f), H* för vinkel 88,90 och 92
H_kritisk = [];
alpha_list = [];
v0_list = [];
figure(4)
% Loop för olika vinklar och plottar
i = 1;
alpha_legend = {};
alpha_sweep = [];
vinklar = [92:-2:88];
for alpha = vinklar
    H_alpha = Bisection(@Diff_Alpha, start, h);                  
    H_kritisk = [H_kritisk; H_alpha];                             
    alpha_list = [alpha_list alpha];                                   
    route_critical_alpha=RK4(h,H_alpha);                           
    trajectory_sweep_alpha=Paths(route_critical_alpha);           
    v0_list = [v0_list trajectory_sweep_alpha.v_sweep];   
    polar(trajectory_sweep_alpha.phi, trajectory_sweep_alpha.r)
    view([90 -90])                      
    hold on
    legend_text = sprintf('%c=%d', 'a', alpha);
    v0_fel(i) = trajectory_sweep_alpha.v_err;
    H_fel(i) = trajectory_sweep_alpha.r_err;
    alpha_legend{i} = legend_text;
    v0(i) = trajectory_sweep_alpha.v_sweep;
    i = i + 1;
end

% Transponat på listor för att göra engöra en tabell.
vinklar = vinklar';
v0 = v0';
v0_fel = v0_fel';
H_fel = H_fel';

% Felberäkningar

title('Banor för olika vinklar')
legend(alpha_legend)

alpha_table = table(vinklar, H_kritisk, H_fel, v0, v0_fel);
disp('Parametrar för de olika vinklarna')
disp(alpha_table)

% Undersöka konvergensen för RK4
disp('Undersöka konvergensen för RK4')
alpha = 90;
error_vector=[];
no_of_tests = 3;

% Loop för olika steglängder
for a=0:no_of_tests-1
    h_test = h*2^a;
    % Utföra RK4 för att sedan göra konvergensberäkning
    u = [H_critical+1 0 0 0];                          % B.V
    t = 0;                                          % start tid
    trajectory_error = struct(    't',        t,...
                        'r',        u(1),...
                        'rd',     u(2),...
                        'phi',      u(3),...
                        'phid',   u(4));
    times = 0:h_test:2;

    for i = 1:length(times)-1
        % RK4
        k1 = h_test*Function(times(i), u, H_critical);
        k2 = h_test*Function(times(i)+0.5*h_test, u+0.5*k1, H_critical);
        k3 = h_test*Function(times(i)+0.5*h_test, u+0.5*k2, H_critical);
        k4 = h_test*Function(times(i)+h_test, u+k3, H_critical);

        u = u + (k1+2*k2+2*k3+k4)/6;
        t = times(i) + h_test; 
        %Stoppa alla värden i en lista
        trajectory_error.t = [trajectory_error.t; t];
        trajectory_error.r = [trajectory_error.r; u(1)];
        trajectory_error.rd = [trajectory_error.rd; u(2)];
        trajectory_error.phi = [trajectory_error.phi; u(3)];
        trajectory_error.phid = [trajectory_error.phid; u(4)];  
    end
    routes_err(a+1) = trajectory_error;
    % Spara alla felberäkningar
    end_value=[ trajectory_error.r(end),...
                trajectory_error.rd(end),...
                trajectory_error.phi(end),...
                trajectory_error.phid(end)];
   error_vector=[error_vector; end_value];
end

r = Table_Data(error_vector(:,1));
rd = Table_Data(error_vector(:,2));
phi = Table_Data(error_vector(:,3));
phid = Table_Data(error_vector(:,4));

tab = table(r, rd, phi, phid);
rows = {'(h)' '(2h)' '(4h)' 'diff(h; 2h)e6' 'diff(2h; 4h)e6' 'quotient' 'relativ fel'};
tab.Properties.RowNames = (rows);
disp(tab)

%Feluppskattning av polynomanpassningen
h = 0.01;
u = [H_critical+1 0 0 0];
trajectory_small = RK4(h, H_critical+1);
trajectory_big = struct(  't',    trajectory_small.t(1:2:end),...
                    'r',    trajectory_small.r(1:2:end),...
                    'rd', trajectory_small.rd(1:2:end),...
                    'phi',  trajectory_small.phi(1:2:end),...
                    'phid',   trajectory_small.phid(1:2:end));

x1 = [];
y1 = [];
x2 = [];
y2 = [];
phi1=[];
phi2=[];

k = 2;
%INterpolation för 2h
for n = 1:length(trajectory_big.t) - 1 
    [x, y, phi] = Hermite(trajectory_big, n, k);
    x1 = [x1 x];
    y1 = [y1 y];
    phi1=[phi1 phi];
    hold on
end

k = 1;
%interpolation gör h
for n = 1:length(trajectory_small.t) - 1 
    [x, y, phi] = Hermite(trajectory_small, n, k);
    x2 = [x2 x];
    y2 = [y2 y];
    phi2=[phi2 phi];
    hold on
end
%Absoluta fel
err_abs = max(abs(y1 - y2));
ind = find(max(abs(y1 - y2)) == err_abs);
err_rel = err_abs/y1(ind);

%Svaret
fprintf('Hermitefel %fe-6\n', err_rel*1e6)
  
function [r_out]=Diff_Alpha(h, H)
    %Funktion för beräkning av kurvorna för olika vinklar
    %Genom RK4

    trajectory=RK4(h,H);

    trajectory_sweep=Paths(trajectory);
    r_out=trajectory_sweep.r_sweep-1; 
end

function [v_vec] = Table_Data(error_list)
e = abs(diff(error_list));
quotient = e(2)/e(1);
rel_err = (error_list(1)-error_list(2))/error_list(1);
%Vektor med reslutaten av konvergensberäkningen
v_vec = error_list;
v_vec(4:5) = e*1e6;
v_vec(6) = quotient;
v_vec(7) = abs(rel_err);
end

function [Lenght] = Path_Lenght(x, C,y)
%  Funktion för att beräkna läångden på kurvan.
%  Formel som används: sqr(1+f'^2), h=>0.
vecL=[];

for i=1:length(x)-1;
    dx=x(i+1)-x(i);
    % Om det är ett polynom, hitta punkterna:
    if y==0
    dy=polyval(C,x(i+1))-polyval(C,x(i));
    % Annars linjär summa. 
    elseif C==0
        dy=y(i+1)-y(i);
    end
dS=sqrt(dx^2+dy^2);
vecL=[vecL dS];
end
Lenght=sum(vecL);
end 
