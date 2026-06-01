% far field di un'onda piana che rimbalza contro una circonferenza 
% calcoli fatti su un foglio di brutta, ho fatto la foto
clear; clc; close all; 

R = 1; % raggio circonferenza
z0 = 0 + 0*1i; % centro circonferenza
k = 20; 
L_max = ceil(k*R)+10;
theta = pi/4; 
d = exp(1i*theta);

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
mask = (X.^2 + Y.^2 < R); % complementare del disco  di raggio R + 0.25

u_scat = zeros(size(Z));

fprintf("Computing u_ref \t")
% calcolo il campo scatterato (non strettamente necessario) e il far field
% pattern (vedere foto sul telefono per la formula usata)
u_coeff = zeros(2*L_max+1, 1);
theta = linspace(0, 2*pi, 360)';
farField = zeros(size(theta));

for l = -L_max:L_max
    u_scat = u_scat - (1i)^l * besselj(l, k*R) ./ besselh(l, 1, k*R) .* besselh(l, 1, k*abs(Z)) .* (Z./(d*abs(Z))).^l;
    u_coeff(l+L_max+1) = (1i)^l * besselj(l, k*R) ./ besselh(l, 1, k*R) * (1/d)^l;
    farField = farField + sqrt( 2/(pi*k) ) * exp(-pi*1i/4) * u_coeff(l+L_max+1) * exp(-l*pi*1i/2) * exp(1i*l*(theta));
    % se sostituisco nella riga sopra theta + pi al posto di theta e
    % considero l'angolo di incidenza dell'onda come theta + pi ho lo
    % stesso plot in polarplot che se non mettessi pi in entrambi i posti
end

u_scat(mask) = complex(NaN, NaN);
fprintf("done\n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
pcolor(X, Y, abs(u_scat)); shading flat; colorbar; axis square; axis off; colormap("hot")

figure;
% plot(-L_max:L_max, real(u_coeff), LineWidth=2, DisplayName="parte reale")
hold on
grid on
% plot(-L_max:L_max, imag(u_coeff), LineWidth=2, DisplayName="parte immaginaria")
plot(-L_max:L_max, abs(u_coeff), LineWidth=2, DisplayName="valore assoluto")
% title("Valore assoluto dei coefficieti del campo scatterato", Interpreter="latex")
% legend(Location="bestoutside")

% Guardando i plot dei valori assoluti dei coefficienti di u_scat si vede
% che per valori di \ell \in \mathbb{Z} molto più grandi o più piccoli (in
% modulo) del valore k*R si hanno coefficienti molto piccoli (in modulo)

figure; 
minFarField = min(log10(abs(farField)));

polarplot(theta, -min(log10(abs(farField))) + log10(abs(farField)), lineWidth = 2)
disp(minFarField)



