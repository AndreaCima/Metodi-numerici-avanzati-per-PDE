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

loglog(N, Err, 'bo-', LineWidth=2)
grid on
hold on
title("Error", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("Error", Interpreter="latex")

figure; 
semilogx(N, times(:, 3), 'bo-', LineWidth=2)
grid on
hold on
title("Plotting time", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("time (s)", Interpreter="latex")


%% Confronto con MPSpack
clear; clc; close all; 

N = 2.^(4:12); % N = 16, 32, ...., 4096
T = load("MPSpackBenchmarkSquareScatt.mat");
u_ref = T.u;

mask = linspace(-1.5, 1.5, 151);
mask = mask(1:end-1);
[X_mask, Y_mask] = meshgrid(mask);
mask = ~(X_mask.^2 + Y_mask.^2 <= 1); % disco unitario

Err = zeros(length(N), 1);
times = zeros(length(N), 3);
H = zeros(length(N), 1);

for s = 1:length(N)
    fprintf("N = %d ", N(s))
    [u_scat, time, h] = BEM_func(N(s));
    fprintf("\t h = %f\n", h(1))
    Err(s) = norm(u_scat(mask)-u_ref(mask), 2) / norm(u_ref(mask), 2);
    times(s, :) = time;
    H(s) = h(1);
end

p = polyfit(log(N(6:end)), log(Err(6:end)), 1); 
slope = p(1); % ordine di convergenza

f1 = figure; 
loglog(N, Err, 'bo-', LineWidth=2)
grid on
hold on
title("Error")
xlabel("Degrees of freedom")
ylabel("Error")

f2 = figure; 
loglog(N, times(:, 3), 'bo-', LineWidth=2)
grid on
hold on
title("Plotting time", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("time (s)", Interpreter="latex")

if ~exist('figures','dir')
    mkdir('figures');
end

exportgraphics(f1, 'figures/erroreMPSpack.png', 'Resolution', 300);
exportgraphics(f2, 'figures/timeplotMPSpack.png', 'Resolution', 300);








