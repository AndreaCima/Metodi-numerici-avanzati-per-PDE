function [u_scat, u_tot, times, h_k] = BEM_curv(N, obs, obs_der, u_inc, options)

arguments 
    N
    obs
    obs_der
    u_inc
    options.k = 20
    options.x_lim = [-2 2]; 
    options.y_lim = [-2 2];
    options.n_points_plot = 300
end


k = options.k;
x_lim = options.x_lim;
y_lim = options.y_lim;
n_points_plot = options.n_points_plot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETRI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs_der_abs = @(t) abs(obs_der(t)); 

n_gauss_pts_plot = 5;
n_gauss_pts_off_diag = 5;
n_gauss_pts_on_diag = 30;

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GEOMETRIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = linspace(0, 2*pi, N+1)';
L = integral(obs_der_abs, 0, 2*pi); % perimetro
h_k = L/N; % ampiezza media degli elementi
x_k = obs( (t(1:N) + t(2:N+1)) / 2 ); % punti medi

%%%%%%%%%%%%%%%%%%%%%%%%%%  ASSEMBLAGGIO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
F = -u_inc(x_k);

A = zeros(N, N);
[xq_off_diag, wq_off_diag] = gaussquad(n_gauss_pts_off_diag);

dt = 2*pi/N; 
t_q = t(1:N) + dt*xq_off_diag.';
y_q = obs(t_q);
jac_q = obs_der_abs(t_q);

for j = 1:N
    r = abs(x_k(j)-y_q);
    integrand = besselh(0, 1, k*r);

    A(j, :) = ( (1i/4) * (integrand .* jac_q) * wq_off_diag * dt ).';
end

% fix diagonale di A
[xq_on_diag, wq_on_diag] = gaussquad(n_gauss_pts_on_diag);
t_mid = (t(1:N) + t(2:N+1)) / 2;
h_j = dt * obs_der_abs(t_mid);

for j = 1:N
    y_q_diag = (h_j(j)/2) * xq_on_diag; 
    integrand = besselh(0, 1, k*y_q_diag);
    A(j, j) = (1i/2)*(h_j(j)/2)*wq_on_diag.'*integrand;
end

time_assembling = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%  SISTEMA LINEARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
psi = A\F; 
time_lin_sist = toc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
x_plot = linspace(x_lim(1), x_lim(2), n_points_plot+1);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot+1);
[X, Y] = meshgrid(x_plot(1:end-1), y_plot(1:end-1));
Z = X + 1i*Y;

u_scat = zeros(size(Z));
xv = real(obs(t(1:N)));
yv = imag(obs(t(1:N)));
in = inpolygon(X, Y, xv, yv);

[xq_plot, wq_plot] = gaussquad(n_gauss_pts_plot);
t_q_plot = t(1:N) + dt*xq_plot.';
y_q_plot = obs(t_q_plot);
jac_q_plot = obs_der_abs(t_q_plot);

for j = 1:numel(Z)
    if ~in(j)
        r = abs(Z(j)-y_q_plot);
        integrand = besselh(0, 1, k*r);
        u_scat(j) = sum( (1i/4) * integrand .* psi .* jac_q_plot * dt * wq_plot );
    end
end
time_plot = toc; 

times = [time_assembling, time_lin_sist, time_plot];

u_inc_grid = u_inc(Z);
u_tot = u_scat + u_inc_grid;

u_tot(in) = complex(NaN, NaN);
u_inc_grid(in) = complex(NaN, NaN);
u_scat(in) = complex(NaN, NaN);

end