clear; clc; close all; 
k = 5;
x0 = 1 + 1i;


u_circ = @(x) exp(1i * k * abs(x - x0));

x_plot = linspace(real(x0)-5, real(x0)+5, 400);
y_plot = linspace(imag(x0)-5, imag(x0)+5, 400);
[X,Y] = meshgrid(x_plot, y_plot);
Z = X + 1i*Y;

u = u_circ(Z);

figure; 
pcolor(X, Y, real(u)); shading flat
colorbar
title('Onda circolare, parte reale')

figure; 
pcolor(X, Y, imag(u)); shading flat
colorbar
title('Onda circolare, parte immaginaria')

figure; 
pcolor(X, Y, abs(u)); shading flat
colorbar
title('Onda circolare, valore assoluto')





