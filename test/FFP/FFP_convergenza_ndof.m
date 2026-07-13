% studio errore FFP rispetto ad una soluzione di riferimento (calcolata con
% un elevato numero di gradi di libertà) all'aumentare dei gradi di libertà

clear; clc; close all; 
max_n = 13;
N_ref = 2^max_n;
k = 20;
v = 0.5*[1+1i, -1+1i, -1-1i, 1-1i];
theta = pi/4;
fprintf('Computing FFP_ref \t')
FFP_ref = Far_Field(N_ref, k, v, theta);
fprintf('done \n')

N = 2.^(4:max_n-1);
Err = zeros(length(N), 1);
H = zeros(length(N), 1);


for s = 1:length(N)
    fprintf('N = %d \t', N(s))
    [FFP, h] = Far_Field(N(s), k, v, theta);
    fprintf('done\n')
    Err(s) = norm(FFP - FFP_ref, 2) / norm(FFP_ref, 2);
    H(s) = h;
end


figure; 
loglog(1./N, Err, 'bo-', LineWidth=2)
grid on

p = polyfit(log(H), log(Err), 1); 
slope = p(1); % ordine di convergenza
clc
fprintf('Slope = %.4f\n', slope)

