% Studio il comportamento della densità psi al variare dei gradi di libertà e
% confrontando il risultato con una soluzione ottenuta con una mesh molto
% fine
clear; clc; close all; 


% parametri
v = 0.5*[-1+1i, -1-1i, 1-1i, 1+1i]; 
k = 20;
n_max = 12;
N_ref = 2^n_max; 
theta = -pi/4;
N = 2.^(4:(n_max-1));


fprintf("Computing psi_ref \t")
[psi_ref, ~, ~] = BEM_psi(N_ref, k, v, theta);
fprintf("done\n")

H = zeros(length(N), 1);
Err = zeros(length(N), 1);
times = zeros(length(N), 1);

for s = 1:length(N)

    fprintf("N = %i\t", N(s))
    [psi, time, h] = BEM_psi(N(s), k, v, theta);
    psi = repelem(psi, 2^(length(N) - s + 1));

    fprintf("h = %.3f\n", h)

    Err(s) = norm( psi - psi_ref, 2) ./ norm(psi_ref, 2);
    H(s) = h;
    times(s) = time(2);


end


p = polyfit(log10(H), log10(Err), 1);
slope = p(1);
fprintf("Slope = %f\n", slope)

figure; 
loglog(N, times, 'bo-', LineWidth=2)
grid on; 
title("Assembling time vs Ndof", Interpreter="latex")
xlabel("degrees of freedom", Interpreter="latex")
ylabel("time (s)")

figure; 
loglog(N, Err, 'bo-', LineWidth=2)
grid on;
title("Error vs Ndof", Interpreter="latex")
xlabel("degrees of freedom", Interpreter="latex")
ylabel("error")

figure; 
loglog(H, Err, 'bo-', LineWidth=2)
grid on; hold on;
title("Error vs mesh size", Interpreter="latex")
xlabel("mesh size", Interpreter="latex")
ylabel("error")

