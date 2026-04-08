clear; clc; close all; 

% Parameters
v = [1+1i, -1+1i, -1-1i]; % vertices of \Gamma (counterclockwise)
k = 5;
N = 50; % degrees of freedom
x_lim = [-2 2];
y_lim = [-2 2];

% Geometry
v = [v v(1)];
n_sides = length(v)-1;
perimeter = sum(abs(v - [v(2:end) v(1)]));
side_length = abs(diff(v));
side_percent = side_length./perimeter;
assert(sum(side_percent)==1);

N_side = ceil(N*side_percent);
p_k = cell(n_sides, 1); % extreme points of an element
x_k = cell(n_sides, 1); % collocation points (midpoints)
tau_k = zeros(n_sides, 1); % tangent vectors (I compute a single one for each side )

h_k = zeros(n_sides, 1); % mesh size per element, elements on the same side share the same value

for i = 1:n_sides
    p_k{i} = linspace(v(i), v(i+1), N_side(i));
    x_k{i} = ( p_k{i}(1:end-1) + p_k{i}(2:end) ) ./2;
    tau_k(i) = (p_k{i}(2)-p_k{i}(1)) / abs(p_k{i}(2)-p_k{i}(1));
    h_k(i) = abs(p_k{i}(2)-p_k{i}(1));
end

plot(real(v), imag(v), LineWidth=2)
grid on
hold on
for i = 1:n_sides
    scatter(real(p_k{i}), imag(p_k{i}), 20, 'filled', 'ro')
    scatter(real(x_k{i}), imag(x_k{i}), 20, 'filled', 'green', 'd')
end

xlim(x_lim)
ylim(y_lim)