clear; clc; close all; 
% calcolo una soluzione molto accurata con il metodo BEM (N alto) e poi la confronto
% con soluzioni progressivamente meno accurate per vedere come si comporta
% l'errore
N_ref = 200;
N_val = 160:-20:20;
n_test = length(N_val);

Err = zeros(n_test, 1);
times = zeros(n_test, 3);

[u_scat_ref, time_assembling_ref, time_lin_sist_ref, time_plot_ref] = BEM_func(N_ref);
for s=1:n_test
    disp(s)
    [u_scat, time_assembling, time_lin_sist, time_plot] = BEM_func(N_val(s));
    mask = ~isnan(u_scat) & ~isnan(u_scat_ref);
    Err(s) = norm(u_scat(mask)-u_scat_ref(mask), 2);
    times(s, :) = [time_assembling, time_lin_sist, time_plot];
end

plot(N_val, Err, LineWidth=2)
grid on
hold on
title("Error")
xlabel("Degrees of freedom")
ylabel("Error")

