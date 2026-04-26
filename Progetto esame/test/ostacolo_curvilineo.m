% metodo BEM con onda piana (per il momento) che sbatte contro ostacolo curvilineo
clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
z0 = 0+0*1i; % centro della circonferenza
R = 1; % raggio

% obs = @(t) z0 + R*exp(1i*t); % circonferenza unitaria, rappresenta il mio ostacolo
% obs_der = @(t) 1i*R*exp(1i*t); % derivata
% obs_der_abs = @(t) abs(1i*R*exp(1i*t)); % abs(derivata), mi serve per sapere il valore del perimetro

obs = @(t) (cos(t) + 0.65*(cos(2*t) - 1)) + 1i*(1.5*sin(t));
obs_der = @(t) (-sin(t) - 1.3*sin(2*t)) + 1i*(1.5*cos(t));
obs_der_abs = @(t) abs(obs_der(t));


k = 20;
N = 500; % degrees of freedom. (Usually N ~ k ~ 1/h). Here I use N=4*k or N = 1000 for some test
x_lim = [-2 2]; 
y_lim = [-2 2];
theta = pi;
n_gauss_pts_plot = 5;
n_gauss_pts_off_diag = 5;
n_gauss_pts_on_diag = 30;
n_points_plot = 300;
u_inc=@(x) exp(1i * k * real(x*exp(-1i*theta))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry
t = linspace(0, 2*pi, N+1)';
L = integral(obs_der_abs, 0, 2*pi); % perimetro
h_k = L/N; % ampiezza elementi
p_k = obs(t);
p_k = p_k(1:N);
x_k = obs( (t(1:N) + t(2:N+1)) / 2 ); % punti medi
tau_k = obs_der(t(1:N));
tau_k = tau_k./abs(tau_k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemblaggio A e F
tic
F = -u_inc(x_k);

A = zeros(N, N);
[xq_off_diag, wq_off_diag] = gaussquad(n_gauss_pts_off_diag);

dt = 2*pi/N; 

y_q = obs(t(1:N) + dt*xq_off_diag.');

for j = 1:N
    r = abs(x_k(j)-y_q);
    integrand = besselh(0, 1, k*r);
    A(j, :) = sum((1i/4) * integrand .* h_k * wq_off_diag, 2);
end

% fix the diagonal of A

[xq_on_diag, wq_on_diag] = gaussquad(n_gauss_pts_on_diag);
y_q = h_k/2 * xq_on_diag;

for j = 1:N
    integrand = besselh(0, 1, k*y_q);
    A(j, j) = (1i/2)*(h_k/2)*wq_on_diag.'*integrand;
end

time_assembling = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sistema lineare
tic
psi = A\F; 
time_lin_sist = toc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
tic

x_plot = linspace(x_lim(1), x_lim(2), n_points_plot+1);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot+1);
[X, Y] = meshgrid(x_plot(1:end-1), y_plot(1:end-1));
Z = X + 1i*Y;

u_scat = zeros(size(Z));
xv = real(obs(t(1:N)));
yv = imag(obs(t(1:N)));

in = inpolygon(X, Y, xv, yv);

[xq_plot, wq_plot] = gaussquad(n_gauss_pts_plot);

y_q = obs(t(1:N) + dt*xq_plot.');

for j = 1:numel(Z)
    if ~in(j)
        r = abs(Z(j)-y_q);
        integrand = besselh(0, 1, k*r);
        u_scat(j) = sum(  (1i/4) * integrand .* psi * h_k * wq_plot);
    end
end

time_plot = toc; 

u_inc_grid = u_inc(Z);
u_tot = u_scat + u_inc_grid;

u_tot(in) = complex(NaN, NaN);
u_inc_grid(in) = complex(NaN, NaN);
u_scat(in) = complex(NaN, NaN);

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

fprintf("Time to assemble A and F = %f seconds\n", time_assembling)
fprintf("Time to solve the linear system = %f seconds\n", time_lin_sist)
fprintf("Time to plot the solution via representation formula = %f seconds\n", time_plot)

function [x, w] = gaussquad(q)
% quadrature nodes and weights for the Gauss quadrature on [0 1]
B = ( 1:(q-1) )./ sqrt( 4*( 1:(q-1) ).^2 -1 );
[V, D] = eig( diag(B, -1) + diag(B, 1) );
x = ( diag(D)+1 )/2;
w = ( V(1, :).*V(1, :) )';

end







