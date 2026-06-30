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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcolo il far field di riferimento

l_target = 0.7;
alpha_target = pi/18;

v = get_vertices(l_target, alpha_target);

[FFP_target, ~] = Far_Field(N, k, v, theta);

clear v

func = @(x) compute_err(x(1), x(2), FFP_target, N, k, theta); % funzione che devo minimizzare
% bounds 
lb = [0.1, 0];
ub = [5, pi/2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ottimizzazione surrogateopt + fminsearch

options_surf = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', 60, ... 
    'Display', 'iter');
% trovo il guess iniziale che poi passo a fminseaech
[x_global, fval_global] = surrogateopt(func, lb, ub, options_surf); 
 
options = optimset('Display', 'iter', 'TolFun', 1e-10);

[x_opt_fmin, err_opt_fmin] = fminsearch(func, x_global, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle swarm

options_swarm = optimoptions('particleswarm', ...
    'SwarmSize', 100, ...            
    'FunctionTolerance', 1e-6, ...  
    'ObjectiveLimit', 1e-8, ...
    'MaxStallIterations', 100, ...   
    'MaxIterations', 2000, ...       
    'Display', 'iter');
% devo specificare il numero di variabili, in questo caso 2 
[x_opt_swarm, err_opt_swarm] = particleswarm(func, 2, lb, ub, options_swarm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stampa dei risultati

metriche = {'Lato'; 'Angolo'; 'Errore'};
target = [l_target; alpha_target; 0]; % il target dell'errore sarebbe 0

% risultati fminsearch
res_fmin = [x_opt_fmin(1); x_opt_fmin(2); err_opt_fmin];

% risultati particleswarm
res_swarm = [x_opt_swarm(1); x_opt_swarm(2); err_opt_swarm];


tabella_confronto = table(metriche, target, res_fmin, res_swarm, ...
    'VariableNames', {'Parametro', 'Target', 'Fmin', 'Swarm'});

disp(tabella_confronto)



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