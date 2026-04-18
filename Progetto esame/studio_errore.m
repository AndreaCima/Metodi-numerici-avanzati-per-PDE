clear; clc; close all; 
% calcolo una soluzione molto accurata con il metodo BEM (N alto) e poi la confronto
% con soluzioni progressivamente meno accurate per vedere come si comporta
% l'errore
N_ref = 150;
N_val = 90:-10:20;

Err = zeros(length(N_val), 1);
times = zeros(length(N_val), 3);

[u_scat_ref, ~] = BEM_func(N_ref);

for s = 1:length(N_val)
    fprintf("N = %d\n", N_val(s))
    [u_scat, time] = BEM_func(N_val(s));
    mask = ~isnan(u_scat) & ~isnan(u_scat_ref);
    Err(s) = norm(u_scat(mask)-u_scat_ref(mask), 2);
    times(s, :) = time;
end

semilogy(N_val, Err, LineWidth=2)
grid on
hold on
title("Error")
xlabel("Degrees of freedom")
ylabel("Error")


%% Confronto con MPSpack
clear; clc; close all; 

N = 16*(1:256);
n = load("MPSpackBenchmarkSquareScatt.mat");

mask = linspace(-1.5, 1.5, 200);
[X_mask, Y_mask] = meshgrid(mask);
mask_ball = (X_mask.^2 + Y_mask.^2 <= 1); % disco unitario

for n = N
    fprintf("N = %d", n)
    []
end






