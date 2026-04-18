clear; clc; close all; 
% calcolo una soluzione molto accurata con il metodo BEM (N alto) e poi la confronto
% con soluzioni progressivamente meno accurate per vedere come si comporta
% l'errore
N_ref = 300;
N = 20:10:200;

Err = zeros(length(N), 1);
times = zeros(length(N), 3);
H = zeros(length(N), 1);

[u_scat_ref, ~, ~] = BEM_func(N_ref);

mask = ~isnan(u_scat_ref);

for s = 1:length(N)
    fprintf("N = %d", N(s))
    [u_scat, time, h] = BEM_func(N(s));
    fprintf("\t h = %f\n", h(1))
    H(s) = h(1);

    Err(s) = norm(u_scat(mask)-u_scat_ref(mask), 2);
    times(s, :) = time;
end

loglog(N, Err, LineWidth=2)
grid on
hold on
title("Error", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("Error", Interpreter="latex")

figure; 
semilogx(N, times(:, 3), LineWidth=2)
grid on
hold on
title("Plotting time", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("time (s)", Interpreter="latex")


%% Confronto con MPSpack
clear; clc; close all; 

N = 16*(1:40);
T = load("MPSpackBenchmarkSquareScatt.mat");

mask = linspace(-1.5, 1.5, 150);
[X_mask, Y_mask] = meshgrid(mask);
mask = ~(X_mask.^2 + Y_mask.^2 <= 1); % disco unitario

Err = zeros(length(N), 1);
times = zeros(length(N), 3);
H = zeros(length(N), 1);

for s = 1:length(N)
    fprintf("N = %d ", N(s))
    [u_scat, time, h] = BEM_func(N(s));
    fprintf("\t h = %f\n", h(1))
    Err(s) = norm(u_scat(mask)-T.ui(mask), 2) / norm(T.ui(mask));
    times(s, :) = time;
    H(s) = h(1);
end
figure; 
loglog(N, Err, LineWidth=2)
grid on
hold on
title("Error")
xlabel("Degrees of freedom")
ylabel("Error")







