%Confronto con MPSpack
clear; clc; close all; 

N = 2.^(4:12); % N = 16, 32, ...., 4096
T = load("MPSpackBenchmarkSquareScatt.mat");
u_ref = T.u;
k = 20;
theta = -pi/4;
v = 0.5* [1+1i, -1+1i, -1-1i, 1-1i];

mask = linspace(-1.5, 1.5, 151);
mask = mask(1:end-1);
[X_mask, Y_mask] = meshgrid(mask);
mask = ~(X_mask.^2 + Y_mask.^2 <= 1); % complementare del disco unitario

Err = zeros(length(N), 1);
times = zeros(length(N), 3);
H = zeros(length(N), 1);

for s = 1:length(N)
    fprintf("N = %d ", N(s))
    [u_scat, time, h] = BEM_func(N(s), k, v, theta);
    fprintf("\t h = %f\n", h)
    Err(s) = norm(u_scat(mask)-u_ref(mask), 2) / norm(u_ref(mask), 2);
    times(s, :) = time;
    H(s) = h;
end

p = polyfit(log(H(end-3:end)), log(Err(end-3:end)), 1); 
slope = p(1); % ordine di convergenza

f1 = figure; 
loglog(N, Err, 'bo-', LineWidth=2)
grid on
hold on
title("Error vs Ndof", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("Error", Interpreter="latex")

f2=figure; 
loglog(H, Err, 'bo-', LineWidth=2, DisplayName='Collocazione')
hold on
loglog(logspace(-3, -2, 200), logspace(-3, -2, 200).^(4/3), 'k--', 'LineWidth', 1.5, DisplayName='$h^{4/3}$')
legend(Location='northwest', Interpreter='latex')
grid on 


title("Error vs mesh size", Interpreter="latex")
xlabel("Mesh size", Interpreter="latex")
ylabel("Error", Interpreter="latex")


f3 = figure; 
loglog(N, times(:, 3), 'bo-', LineWidth=2)
grid on
hold on
title("Plotting time", Interpreter="latex")
xlabel("Degrees of freedom", Interpreter="latex")
ylabel("time (s)", Interpreter="latex")

dirpath = 'figures';
if (exist(dirpath, 'dir') == 0), mkdir(dirpath); end

% exportgraphics(f1, 'figures/Err_N_MPSpack.png', 'Resolution', 300);
% exportgraphics(f2, 'figures/Err_mesh_MPSpack.png', 'Resolution', 300);
% exportgraphics(f3, 'figures/timeplotMPSpack.png', 'Resolution', 300);








