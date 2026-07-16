clear; clc; close all;

v = [ ...
    -3+0i,  3+0i,  3+1i,  1.8+1i, ...
    1.8+6i, 2.8+6i, 2.8+8i, 1.6+8i, ...
    1.6+6.8i, 0.6+6.8i, 0.6+8i, -0.6+8i, ...
    -0.6+6.8i, -1.6+6.8i, -1.6+8i, -2.8+8i, ...
    -2.8+6i, -1.8+6i, -1.8+1i, -3+1i ...
];

k = 20;
N = 1000; % degrees of freedom. (Usually N ~ k ~ 1/h). Here I use N=4*k or N = 1000 for some test
x_lim = [-5 5]; 
y_lim = [0 10];
theta = pi/2;
n_gauss_pts_plot = 5;
n_gauss_pts_off_diag = 5;
n_gauss_pts_on_diag = 30;
n_points_plot = 600;
u_inc=@(x) exp(1i * k * real(x*exp(-1i*theta))); 

% Geometry
v = [v v(1)];
n_sides = length(v)-1;
perimeter = sum(abs(v - [v(2:end) v(1)]));
side_length = abs(diff(v));
side_percent = side_length./perimeter;
assert(sum(side_percent) - 1 < 1e-4);

N_side = ceil(N*side_percent);
N = sum(N_side); % new N
p_k = cell(n_sides, 1); % extreme points of an element

x_k = cell(n_sides, 1); % collocation points (midpoints)
tau_k = cell(n_sides, 1); % tangent vectors (I compute a single one for each side )

h_k = cell(n_sides, 1); % mesh size per element, elements on the same side share the same value

for i = 1:n_sides
    p_k{i} = linspace(v(i), v(i+1), N_side(i)+1);
    x_k{i} = ( p_k{i}(1:end-1) + p_k{i}(2:end) ) ./2;
    tau_k{i} = (p_k{i}(2)-p_k{i}(1)) / abs(p_k{i}(2)-p_k{i}(1))*ones(1, length(p_k{i})-1);
    h_k{i} = abs(p_k{i}(2)-p_k{i}(1))*ones(1, length(p_k{i})-1);
end
x_k = [x_k{:}].'; 
p_k = [p_k{:}].'; 
% delete the edges of Gamma (where diff(p_k)=0. I mantain th first one st p_0 = p_N
p_k = p_k(diff(p_k) ~= 0); 
h_k = [h_k{:}].';
tau_k = [tau_k{:}].';

% Assembling A and F
tic
F = -u_inc(x_k);

A = zeros(N, N);
[xq_off_diag, wq_off_diag] = gaussquad(n_gauss_pts_off_diag);
y_q = p_k + h_k * xq_off_diag.' .* tau_k;
for j = 1:N
    
    r = abs(x_k(j)-y_q);
    integrand = besselh(0, 1, k*r);
    A(j, :) = sum((1i/4) * integrand .* h_k * wq_off_diag, 2);
end
% fix the diagonal of A
[xq_on_diag, wq_on_diag] = gaussquad(n_gauss_pts_on_diag);
for j = 1:N
    y_q = h_k(j) * xq_on_diag;
   
    integrand = besselh(0, 1, k*y_q);
    A(j, j) = (1i/2)*(h_k(j)/2)*wq_on_diag.'*integrand;
end
time_assembling = toc;

% Solve the linear system
tic
psi = A\F;
time_lin_sist = toc; 

% Plot the solution
tic
x_plot = linspace(x_lim(1), x_lim(2), n_points_plot+1);
y_plot = linspace(y_lim(1), y_lim(2), n_points_plot+1);
[X, Y] = meshgrid(x_plot(1:end-1), y_plot(1:end-1));
Z = X + 1i*Y;

u_scat = zeros(size(Z));
in = inpolygon(X, Y, real(v), imag(v));
[xq_plot, wq_plot] = gaussquad(n_gauss_pts_plot);

for j = 1:numel(Z)
    if ~in(j)
        y_q = p_k + (h_k * xq_plot.') .* tau_k;
        r = abs(Z(j)-y_q);
        integrand = besselh(0, 1, k*r);
        u_scat(j) = sum(  (1i/4) * integrand .* psi .* h_k * wq_plot);
    end
end
time_plot = toc;

u_inc_grid = u_inc(Z);
u_tot = u_scat + u_inc_grid;

u_tot(in) = complex(NaN, NaN);
u_inc_grid(in) = complex(NaN, NaN);
u_scat(in) = complex(NaN, NaN);

% u_inc
figure; 
pcolor(X,Y,real(u_inc_grid)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{R}u_{inc}$',Interpreter='latex')

figure; 
pcolor(X,Y,imag(u_inc_grid)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{I}u_{inc}$',Interpreter='latex')

figure; 
pcolor(X,Y,abs(u_inc_grid)); shading flat; axis equal; colormap(hot); axis off; colorbar
title('$|u_{inc}|$',Interpreter='latex')

% u_scat
figure; 
pcolor(X,Y,real(u_scat)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{R}u_{scat}$',Interpreter='latex')

figure; 
pcolor(X,Y,imag(u_scat)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{I}u_{scat}$', Interpreter='latex')

figure; 
pcolor(X,Y,abs(u_scat)); shading flat; colormap(hot); axis equal; axis off; colorbar
title('$|u_{scat}|$',Interpreter='latex')

% u_tot
figure; 
pcolor(X,Y,real(u_tot)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{R}u_{tot}$',Interpreter='latex')

figure; 
pcolor(X,Y,imag(u_tot)); shading flat; axis equal; axis off; colorbar
title('$\mathcal{I}u_{tot}$',Interpreter='latex')

figure; 
pcolor(X,Y,abs(u_tot)); shading flat; axis equal; colormap(hot); axis off; colorbar
title('$|u_{tot}|$',Interpreter='latex')

fprintf("Time to assemble A and F = %f seconds\n", time_assembling)
fprintf("Time to solve the linear system = %f seconds\n", time_lin_sist)
fprintf("Time to plot the solution via representation formula = %f seconds\n", time_plot)





