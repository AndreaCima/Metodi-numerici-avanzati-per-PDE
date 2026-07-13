% Studio il comportamento dell'errore al variare dei gradi di libertà e
% confrontando il risultato con una soluzione ottenuta con una mesh molto
% fine
clear; clc; close all; 

v = 0.5* [1+1i, -1+1i, -1-1i, 1-1i];
N_ref = 2^12;
N = 2.^(4:11);
k = 20; 
theta = -pi/4;

x_plot = linspace(-1.5, 1.5, 151); % n_points_plot = 150
[X, Y] = meshgrid(x_plot(1:end-1));

mask = ~inpolygon(X, Y, real(v), imag(v)); % caso generale 

% toglere commento alle mask sotto per confronto con MPSpack
% mask = linspace(-1.5, 1.5, 151);
% mask = mask(1:end-1);
% [X_mask, Y_mask] = meshgrid(mask);
% mask = ~(X_mask.^2 + Y_mask.^2 <= 1); % complementare del disco unitario

fprintf("Computing u_ref \t")
[u_ref, ~, ~] = BEM_func(N_ref, k, v, theta);
fprintf("done \n")

H = zeros(length(N), 1);
Err = zeros(length(N), 1);
times = zeros(length(N), 1);

for s = 1:length(N)
    fprintf("N = %i \t", N(s))
    [u_scat, time, h] = BEM_func(N(s), k, v, theta);
    fprintf("h = %.3f\n", h)

    H(s) = h;
    times(s) = time(3);
    Err(s) = norm(u_scat(mask) - u_ref(mask), 2)./norm(u_ref(mask), 2);

end
p = polyfit(log10(H(end-3:end)), log10(Err(end-3:end)), 1);
slope = p(1);
clc; 
fprintf("Slope = %f\n", slope)


figure; 
loglog(N, times, 'bo-', LineWidth=2)
grid on; 
title("Plotting time vs Ndof", Interpreter="latex")
xlabel("degrees of freedom")
ylabel('time (s)')

figure; 
loglog(N, Err, 'bo-', LineWidth=2)
grid on;
title("Error vs Ndof", Interpreter="latex")
xlabel("degrees of freedom")

figure; 
loglog(H, Err, 'bo-', LineWidth=2)
grid on; hold on;
title("Error vs mesh size", Interpreter="latex")
xlabel("mesh size")



