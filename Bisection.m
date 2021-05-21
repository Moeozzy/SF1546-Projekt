function [root] = Bisection (fhandle, start, h)

h_min = 0.0001;
iterations = 100;
x0 = start(1); 
x1 = start(2); 
it = 1;       % Starting values, it = iteration counter
[eval] = feval(fhandle, h, x0); % Evaluate function at starting value 1
f0 = eval;
h_max = 1;

% Bisection method
while abs(h_max) > h_min
    % Övre gräns
    if it > iterations
        disp('uppnåt övregräns')
        break
    end
    [eval] = feval(fhandle, h, x1);            % Evaluera funktionen
    f1=eval;
    h_max = (x1 - x0)/(f1 - f0)*f1;        
    x0 = x1; 
    f0 = f1;                   % Ny gissning
    x1 = x1 - h_max;                       
    it = it + 1;                       
end

root = x1;