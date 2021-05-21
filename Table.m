function [] = Table (routes)
%Funktion för att skapa en tabell för att det ska vara mer tydligt.
amount_of_rows = 4;

%Fyller tabellen
for i = 1:amount_of_rows
    route = routes(i);
    H(i, 1) = route.H;
    t_sweep(i, 1) = route.t_sweep; 
    t_err(i, 1) = route.t_err;
    r_sweep(i, 1) = route.r_sweep;
    r_err(i, 1) = route.r_err;
    phi_sweep(i, 1) = route.phi_sweep;
    phi_err(i, 1) = route.phi_err;
    v_sweep(i, 1) = route.v_sweep;
    v_err(i, 1) = route.v_err;
end

T = table(H, t_sweep, t_err, r_sweep, r_err, phi_sweep, phi_err);
disp('Värden på parametrarna för raketen')
disp(T)