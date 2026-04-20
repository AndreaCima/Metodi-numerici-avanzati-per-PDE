clear; clc; close all; 
k = 5;
x0 = 1 + 1i;

u_circ = @(x) besselh(0,1, k*abs(x - x0));

x_plot = linspace(real(x0)-2, real(x0)+2, 400);
y_plot = linspace(imag(x0)-2, imag(x0)+2, 400);
[X,Y] = meshgrid(x_plot, y_plot);
Z = X + 1i*Y;

u = u_circ(Z);

figure; 
pcolor(X, Y, real(u)); shading flat
colorbar
title('Soluzione fondamentale, parte reale', Interpreter='latex')

figure; 
pcolor(X, Y, imag(u)); shading flat
colorbar
title('Soluzione fondamentale, parte immaginaria', Interpreter='latex')

figure; 
pcolor(X, Y, abs(u)); shading flat
colorbar
title('Soluzione fondamentale, valore assoluto', Interpreter='latex')





