clear; clc; close all; 

% parameters
v = 0.5*[-1+1i, -1-1i, 1-1i]; 
k = 20;
N = 4*k; 
x0 = .5 + .5i;
x_lim = [real(x0)-2,  real(x0)+2]; 
y_lim = [imag(x0)-2, imag(x0)+2];

n_gauss_pts_plot = 5;
n_gauss_pts_off_diag = 4;
n_gauss_pts_on_diag = 10;
n_points_plot = 200;

% due possibili input
l = 3;
u_inc = @(x) besselh(l,1,k*abs(x-x0)) .* exp(1i*l*angle(x-x0)); 
% u_inc = @(x) (1i/4)*besselh(0,1, k*abs(x - x0)); % sol fondamentale 

[u_scat, u_tot, times, mean_h] = BEM_func(N, k, v, u_inc);

x_plot = linspace(x_lim(1), x_lim(2), 150);
y_plot = linspace(y_lim(1), y_lim(2), 150);
[X, Y] = meshgrid(x_plot, y_plot);
Z = X + 1i*Y;
u_inc_grid = u_inc(Z);

in = inpolygon(X, Y, real(v), imag(v));

u_tot(abs(u_tot) > 2) = complex(NaN, NaN);

u_inc_grid(in) = complex(NaN, NaN);
u_inc_grid(abs(u_inc_grid) > 2) = complex(NaN, NaN);

% u_inc
figure; 
pcolor(X,Y,real(u_inc_grid)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{R}u_{inc}$',Interpreter='latex')

figure; 
pcolor(X,Y,imag(u_inc_grid)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{I}u_{inc}$',Interpreter='latex')

figure; 
pcolor(X,Y,abs(u_inc_grid)); shading flat; axis equal; colormap(hot); axis off; colorbar
title('$|u_{inc}|$',Interpreter='latex')

% u_scat
figure; 
pcolor(X,Y,real(u_scat)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{R}u_{scat}$',Interpreter='latex')

figure; 
pcolor(X,Y,imag(u_scat)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{I}u_{scat}$', Interpreter='latex')

figure; 
pcolor(X,Y,abs(u_scat)); shading flat; colormap(hot); axis equal; axis off; colorbar
title('$|u_{scat}|$',Interpreter='latex')

% u_tot
figure; 
pcolor(X,Y,real(u_tot)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{R}u_{tot}$',Interpreter='latex')

figure; 
pcolor(X,Y,imag(u_tot)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{I}u_{tot}$',Interpreter='latex')

figure; 
pcolor(X,Y,abs(u_tot)); shading flat; axis equal; colormap(hot); axis off; colorbar
title('$|u_{tot}|$',Interpreter='latex')




