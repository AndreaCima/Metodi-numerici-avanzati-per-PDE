clear; clc; close all;

% Considero il far field di una circonferenza di raggio R e guardo il far
% field di poligoni regolari inscritti con numero di lati variabile 
% mi aspetto che aumentando il numero di lati abbia approssimazione sempre più
% accurata

max_edges = 12; 
edges = 2.^(4:max_edges); 

R = 1; % raggio circonferenza
z0 = 0 + 0*1i; % centro circonferenza
k_val = [5, 10, 15, 20]; 
theta_inc = pi/4; 

d = exp(1i*theta_inc);
theta = linspace(0, 2*pi, 360)';
FFP_ref = cell(length(k_val), 1);
for k = 1:length(k_val)
    L_max = ceil(k*R)+20;
    u_coeff = zeros(2*L_max+1, 1);
    FFP_ref{k} = zeros(length(theta), 1);
    for l = -L_max:L_max
        u_coeff(l+L_max+1) = -(1i)^l * besselj(l, k*R) ./ besselh(l, 1, k*R) * (1/d)^l; 
        FFP_ref{k} = FFP_ref{k} + sqrt( 2/(pi*k) ) * exp(-pi*1i/4) * u_coeff(l+L_max+1) * exp(-l*pi*1i/2) * exp(1i*l*(theta));
    end
end


Err = zeros(length(edges), length(k_val));

for k = 1:length(k_val)
    fprintf("Testing con k = %i\n", k_val(k))
    L_max = ceil(k*R)+20;
    for n = progress(1:length(edges))
        v = z0 + R * exp( 1i * (0:edges(n)-1) * 2*pi / edges(n) ); % vertici del poligono inscritto
        [FFP, ~] = Far_Field(edges(n), k, v, theta_inc);
        Err(n, k) = norm(FFP-FFP_ref{k}, 2)/norm(FFP_ref{k}, 2);
    end

end

figure; 
loglog(edges, Err, Marker="o", LineWidth=2)
grid on
hold on 
labels = arrayfun(@(k) sprintf("k = %i", k), k_val, UniformOutput=false);
legend(labels)
xlabel('numero lati', Interpreter = 'latex')
ylabel('errore', Interpreter = 'latex')

