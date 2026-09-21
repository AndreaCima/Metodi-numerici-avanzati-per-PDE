% metodo BEM con onda piana che urta contro ostacolo curvilineo
% plot di parte reale, immaginaria e modulo per campo incidente, scatterato
% e totale 

clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
z0 = 0+0*1i; 
R = 1; 

% obs = @(t) z0 + R*exp(1i*t); % circonferenza unitaria
% obs_der = @(t) 1i*R*exp(1i*t); % derivata
% obs_der_abs = @(t) abs(1i*R*exp(1i*t)); % abs(derivata), mi serve per sapere il valore del perimetro

% kite
obs = @(t) (cos(t) + 0.65*(cos(2*t) - 1)) + 1i*(1.5*sin(t));
obs_der = @(t) (-sin(t) - 1.3*sin(2*t)) + 1i*(1.5*cos(t));
obs_der_abs = @(t) abs(obs_der(t));


k = 20;
N = 1000; 
x_lim = [-2 2]; 
y_lim = [-2 2];
theta = pi;
n_points_plot = 200;

u_inc=@(x) exp(1i * k * real(x*exp(-1i*theta))); 


[u_scat, u_tot, times, h_k] = BEM_curv(N, obs, obs_der, u_inc, 'n_points_plot', n_points_plot);

x_plot = linspace(x_lim(1), x_lim(2), n_points_plot0);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot);
[X, Y] = meshgrid(x_plot, y_plot);
Z = X + 1i*Y;
u_inc_grid = u_inc(Z);

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