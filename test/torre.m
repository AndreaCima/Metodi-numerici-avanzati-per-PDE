clear; clc; close all;

v = [ ...
    -3+0i,  3+0i,  3+1i,  1.8+1i, ...
    1.8+6i, 2.8+6i, 2.8+8i, 1.6+8i, ...
    1.6+6.8i, 0.6+6.8i, 0.6+8i, -0.6+8i, ...
    -0.6+6.8i, -1.6+6.8i, -1.6+8i, -2.8+8i, ...
    -2.8+6i, -1.8+6i, -1.8+1i, -3+1i ...
];

k = 20;
N = 1000; 
x_lim = [-5 5]; 
y_lim = [0 10];
theta = pi/2;
n_gauss_pts_plot = 5;
n_gauss_pts_off_diag = 5;
n_gauss_pts_on_diag = 30;
n_points_plot = 600;
u_inc=@(x) exp(1i * k * real(x*exp(-1i*theta))); 

[u_scat, u_tot, times, mean_h] = BEM_func(N, k, v, u_inc, n_points_plot, x_lim, y_lim);

x_plot = linspace(x_lim(1), x_lim(2), n_points_plot);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot);
[X, Y] = meshgrid(x_plot, y_plot);
Z = X + 1i*Y;
u_inc_grid = u_inc(Z);

in = inpolygon(X, Y, real(v), imag(v));

u_inc_grid(in) = complex(NaN, NaN);


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







