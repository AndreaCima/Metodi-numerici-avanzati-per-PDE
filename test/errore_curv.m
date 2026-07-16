% Errore BEM nel caso in cui l'stacolo sia una circonferenza

clear; clc; close all; 


R = 1; % raggio 
z0 = 0 + 0*1i; % centro

k = 20; 
L_max = ceil(k*R)+10;
theta = pi; 
d = exp(1i*theta);

% parametrizzazione bordo e derivata
obs = @(t) z0 + R*exp(1i*t); 
obs_der = @(t) 1i*R*exp(1i*t); 

x_lim = [-2 2]; 
y_lim = [-2 2];
n_points_plot = 300; 
x_plot = linspace(x_lim(1), x_lim(2), n_points_plot+1);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot+1);
[X, Y] = meshgrid(x_plot(1:end-1), y_plot(1:end-1));
Z = X + 1i*Y; 

u_ref = zeros(size(Z));
fprintf("Computing u_ref \t")
for l = -L_max:L_max
    u_ref = u_ref - (1i)^l * besselj(l, k*R) ./ besselh(l, 1, k*R) .* besselh(l, 1, k*abs(Z)) .* (Z./(d*abs(Z))).^l;
end
fprintf("done\n")

N = 2.^(4:12);
mask = linspace(-2, 2, n_points_plot);
% mask = mask(1:end-1);
[X_mask, Y_mask] = meshgrid(mask);
mask = ~(X_mask.^2 + Y_mask.^2 <= R+0.25); % complementare del disco  di raggio R + 0.25


H = zeros(length(N), 1);
Err = zeros(length(N), 1);
times = zeros(length(N), 3);

for s = 1:length(N)

    t = linspace(0, 2*pi, N(s)+1)';
    xv = real(obs(t(1:N(s))));
    yv = imag(obs(t(1:N(s))));
    in = inpolygon(X, Y, xv, yv); 

    fprintf("N = %d ", N(s))
    [u_scat, time, h] = BEM_curv(N(s), obs, obs_der);
    fprintf("\t h = %f\n", h)

    H(s) = h;
    times(s, :) = time; 

    Err(s) = norm(u_scat(mask)-u_ref(mask), 2) / norm(u_ref(mask), 2);



end

p = polyfit(log(H(end-4:end)), log(Err(end-4:end)), 1); 
slope = p(1); % ordine di convergenza

f1 = figure; 
loglog(N, Err, 'bo-', LineWidth=2)
grid on
hold on
title("Error vs Ndof", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("Error", Interpreter="latex")

f2=figure; 
loglog(H, Err, 'bo-', LineWidth=2, DisplayName='Err')
grid on 
hold on
title("Error vs mesh size", Interpreter="latex")
xlabel("Mesh size", Interpreter="latex")
ylabel("Error", Interpreter="latex")

f3 = figure; 
loglog(N, times(:, 3), 'bo-', LineWidth=2)
grid on
hold on
title("Plotting time", Interpreter="latex")
xlabel("degrees of freedom", Interpreter="latex")
ylabel("time (s)")




