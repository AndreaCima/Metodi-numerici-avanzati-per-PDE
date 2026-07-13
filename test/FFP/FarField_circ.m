% far field di un'onda piana contro una circonferenza 
clear; clc; close all; 

R = 1; 
z0 = 0 + 0*1i; 
k = 40; 
L_max = ceil(k*R)+20;
theta_inc = pi/3; 
d = exp(1i*theta_inc);

obs = @(t) z0 + R*exp(1i*t); 
obs_der = @(t) 1i*R*exp(1i*t); 
obs_der_abs = @(t) abs(1i*R*exp(1i*t)); 


x_lim = [real(z0)-R-2, real(z0)+R+2]; 
y_lim = [imag(z0)-R-2, imag(z0)+R+2];
n_points_plot = 300; 
x_plot = linspace(x_lim(1), x_lim(2), n_points_plot+1);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot+1);
[X, Y] = meshgrid(x_plot, y_plot);
Z = X + 1i*Y; 
mask = (X.^2 + Y.^2 < R); % disco di raggio R 

u_scat = zeros(size(Z));

u_coeff = zeros(2*L_max+1, 1);
theta = linspace(0, 2*pi, 360)';
farField = zeros(size(theta));

for l = -L_max:L_max
    % calcolo il campo scatterato solo per i plot, non serve per il far
    % field
    u_scat = u_scat - (1i)^l * besselj(l, k*R) ./ besselh(l, 1, k*R) .* besselh(l, 1, k*abs(Z)) .* (Z./(d*abs(Z))).^l;
    
    u_coeff(l+L_max+1) = -(1i)^l * besselj(l, k*R) ./ besselh(l, 1, k*R) * (1/d)^l; 
    farField = farField + sqrt( 2/(pi*k) ) * exp(-pi*1i/4) * u_coeff(l+L_max+1) * exp(-l*pi*1i/2) * exp(1i*l*(theta));
    % se sostituisco nella riga sopra theta + pi al posto di theta e
    % considero l'angolo di incidenza dell'onda come theta + pi ho lo
    % stesso plot in polarplot che se non mettessi pi in entrambi i posti
end

u_scat(mask) = complex(NaN, NaN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
pcolor(X, Y, abs(u_scat)); shading flat; colorbar; axis square; axis off; colormap("hot")

figure;
name_plot = sprintf("Modulo coefficenti campo scatterato, kR=%.1f", k*R);
semilogy(-L_max:L_max, abs(u_coeff), LineWidth=2)
title(name_plot, Interpreter="latex")
grid on

% si vede che abs(coeff_\ell) \ll 1 se |\ell| \ge kR + 20

figure; 
minFarField = min(log10(abs(farField)));

polarplot(theta, -min(log10(abs(farField))) + log10(abs(farField)), lineWidth = 2)



