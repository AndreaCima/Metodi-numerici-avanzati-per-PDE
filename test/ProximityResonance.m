clear; clc; close all; 

k = 20;
N = 4*k;
v = 0.5*[0+1i, -1-0i, 0-1i];

n_gauss_pts_plot = 5;
n_gauss_pts_off_diag = 4;
n_gauss_pts_on_diag = 10;
n_points_plot = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry
v = [v v(1)];
n_sides = length(v)-1;
perimeter = sum(abs(v - [v(2:end) v(1)]));
side_length = abs(diff(v));
side_percent = side_length./perimeter;
assert(abs(sum(side_percent)-1)<1e-4);

if max(side_percent) - min(side_percent) < 1e-12
    side_percent = ones(size(side_percent)) / n_sides;
end

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

p_k = p_k(diff(p_k) ~= 0); 
h_k = [h_k{:}].';
tau_k = [tau_k{:}].';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 2*pi/k; 

dist_val = [lambda/4, lambda/2, 3*lambda/4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling matrix A

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


figures = cell(3, 2);

for d = 1:length(dist_val)
    dist = dist_val(d);
    x0 = dist + 0*1i;

    x_lim = [real(x0)-1,  real(x0)+1]; 
    y_lim = [imag(x0)-1, imag(x0)+1];
    
    
    u_inc = @(x) besselh(0,1, k*abs(x - x0)); % sol fondamentale
   
    F = -u_inc(x_k);

    % Solve the linear system
    psi = A\F;
    
    % Plot the solution
    x_plot = linspace(x_lim(1), x_lim(2), n_points_plot);
    y_plot = linspace(y_lim(1), y_lim(2), n_points_plot);
    [X, Y] = meshgrid(x_plot, y_plot);
    Z = X + 1i*Y;
    
    u_scat = zeros(size(Z));
    in = inpolygon(X, Y, real(v), imag(v));
    [xq_plot, wq_plot] = gaussquad(n_gauss_pts_plot);
    
    for j = 1:numel(Z)
        if in(j)
            u_scat(j) = complex(NaN, NaN); % inside the polygon
        else
            val = 0;
            for m = 1:N
                y_q = p_k(m) + (h_k(m)/2) * (xq_plot+1) * tau_k(m);
                r = abs(Z(j) - y_q);
                integrand = besselh(0, 1, k * r);
                val = val + (1i / 4) * psi(m) * (h_k(m)/2)* (wq_plot.' * integrand);
            end
            u_scat(j) = val;
        end
    end
    
    u_inc_grid = u_inc(Z);
    u_inc_grid(in) = complex(NaN, NaN);
    u_tot = u_scat + u_inc_grid;

    figures(:, d) = {u_inc_grid, u_scat, u_tot};

end

% plotting the figures


% u_tot
figure; 
pcolor(abs(figures{3, 1})); shading flat; axis equal; colormap("hot"); axis off; colorbar
title('$dist = \frac{\lambda}{4}$', Interpreter='latex')

figure; 
pcolor(abs(figures{3, 2})); shading flat; axis equal; colormap("hot"); axis off; colorbar
title('$dist = \frac{\lambda}{2}$', Interpreter='latex')

figure; 
pcolor(abs(figures{3, 3})); shading flat; axis equal; colormap("hot"); axis off; colorbar
title('$dist = \frac{3\lambda}{4}$', Interpreter='latex')







