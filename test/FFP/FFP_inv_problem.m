% uso la funzione fminsearch per risolvere un problema inverso in cui dato
% il far field cerco di ricostruire la forma dell'ostacolo. 
% suppongo di sapere che
% è un quadrato centrato nell'origine, ma non conosco nè la lungezza del
% lato, nè di quanto è ruotato

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all; 
% Parametri fissi
N = 2^8;
k = 20;
theta = 0;
tol = 1e-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vettori dei target
l_targets = [0.5, 0.7, 1.5, 3.0];
alpha_targets = [0, pi/18, pi/4, pi/3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parametri degli algoritmi di ottimizzazione
options_surf = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', 200, ... 
    'Display', 'off');

options_fmin = optimset('Display', 'off', ...
    'TolFun', tol);

options_swarm = optimoptions('particleswarm', ...
    'SwarmSize', 100, ...            
    'FunctionTolerance', 1e-6, ...  
    'ObjectiveLimit', tol, ...
    'MaxStallIterations', 100, ...   
    'MaxIterations', 2000, ...       
    'Display', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bounds 
lb = [0.1, 0];
ub = [5, pi/2];

matrix_results = zeros(length(l_targets)*length(alpha_targets), 10); 
row_idx = 1;
for i = 1:length(l_targets)
    for j = 1:length(alpha_targets)

        l_target = l_targets(i);
        alpha_target = alpha_targets(j);

        fprintf('Test in corso: l = %.2f, alpha = %.3f rad ... \n', l_target, alpha_target);

        v = get_vertices(l_target, alpha_target);
        [FFP_target, ~] = Far_Field(N, k, v, theta);

        func = @(x) compute_err(x(1), x(2), FFP_target, N, k, theta);

        % surrogate opt + fmin
        tic;
        [x_global, ~] = surrogateopt(func, lb, ub, options_surf);
        [x_opt_fmin, err_opt_fmin] = fminsearch(func, x_global, options_fmin);
        time_fmin = toc;

        % particle swarm
        tic;
        [x_opt_swarm, err_opt_swarm] = particleswarm(func, 2, lb, ub, options_swarm);
        time_swarm = toc;

        row = [l_target, alpha_target, ...
            x_opt_fmin(1), x_opt_fmin(2), err_opt_fmin, time_fmin, ...
            x_opt_swarm(1), x_opt_swarm(2), err_opt_swarm, time_swarm];

        % row = [l_target, alpha_target, x_opt_fmin(1), x_opt_fmin(2), err_opt_fmin, time_fmin];

        matrix_results(row_idx, :) = row;
        row_idx = row_idx + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cols_names = {'Target_Lato', 'Target_Angolo', ...
    'Fmin_Lato', 'Fmin_Angolo', 'Fmin_Errore', 'Fmin_Tempo_s', ...
    'Swarm_Lato', 'Swarm_Angolo', 'Swarm_Errore', 'Swarm_Tempo_s'};

table_results = array2table(matrix_results, 'VariableNames', cols_names);
disp(table_results)

writetable(table_results, 'inv_problem.csv')



function v = get_vertices(l, alpha)
% data la lungezza del lato e l'orientazione, calcola i vertici del
% quadrato 

v = (l/sqrt(2))* exp( 1i*( (0:3)*pi/2 + pi/4 + alpha ) );

end

function err = compute_err(l, alpha, FFP_target, N, k, theta)

v = get_vertices(l, alpha);
[FFP, ~] = Far_Field(N, k, v, theta);

err = norm(FFP-FFP_target, 2) / norm(FFP_target, 2);
end