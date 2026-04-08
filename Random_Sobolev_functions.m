clear; clc; close all; 
rng(100)
ell_max = 100;
theta = linspace(0, 2*pi, 200);
epsilon = 0.1;
s_val = [0 0.5 1];

colors = ['r', 'g', 'b'];

for idx = 1:length(s_val)
    s = s_val(idx);
    ell = 1:ell_max;
    mod_v = ell.^(-s + epsilon);
    phi_pos = 2*pi*rand(1, ell_max);
    phi_neg = 2*pi*rand(1, ell_max);

    v_pos = mod_v .* exp(1i*phi_pos);
    v_neg = mod_v .* exp(1i*phi_neg);
    v_coeff = [v_neg, 0, v_pos];

    ell_full = -ell_max:ell_max;

    v = zeros(size(theta));
    for k = 1:length(ell_full)
        v = v + v_coeff(k) .* exp(1i * ell_full(k) * theta);
    end

    v_real = real(v);
    v_imag = imag(v);

    subplot(1, 2, 1)
    plot(theta, v_real, colors(idx), 'DisplayName', ['s = ' num2str(s)])
    grid on
    hold on
    subplot(1, 2, 2)
    plot(theta, v_imag, colors(idx), 'DisplayName', ['s = ' num2str(s)])
    grid on
    hold on
end

subplot(1, 2, 1)
title('Real part', 'Interpreter', 'latex')
xlabel('Theta')
ylabel('v_{real}')
legend(Location= "best")
axis tight

subplot(1, 2, 2)
title('Imaginary part', 'Interpreter', 'latex')
xlabel('Theta')
ylabel('v_{imag}')
legend(Location= "best")
axis tight




