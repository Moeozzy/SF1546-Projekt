function [y_vec, err]=Linearpol (x,y,xq)
% Linear interpolation

% Beräkna interpolerade punkterna
y_vec = y(2)+(y(3)-y(2))/(x(3)-x(2))*(xq-x(2));

% Samma med dubbla h
y_vec_2h = y(1)+(y(2)-y(1))/(x(2)-x(1))*(xq-x(1));

% Felberäkning
err = abs(y_vec - y_vec_2h);