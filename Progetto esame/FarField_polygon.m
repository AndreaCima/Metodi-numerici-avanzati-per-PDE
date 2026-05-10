% Far field pattern di un'onda piana che rimbalza contro un poligono
clear; clc; close all; 

k = 20; 
N = 8*k;

v = [1+0*1i, 0+1i, 0+0*1i]; % vertici (senso antiorario)

theta = pi/3; 
u_inc=@(x) exp(1i * k * real(x*exp(-1i*theta))); 

n_gauss_pts_farField = 5;
n_gauss_pts_off_diag = 5;
n_gauss_pts_on_diag = 30;

% Geometry
v = [v v(1)];
n_sides = length(v)-1;
perimeter = sum(abs(v - [v(2:end) v(1)]));
side_length = abs(diff(v));
side_percent = side_length./perimeter;
assert(abs( sum(side_percent)-1 ) < 1e-4);

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

% delete the edges of Gamma (where diff(p_k)=0)
p_k = p_k(diff(p_k) ~= 0); 
h_k = [h_k{:}].';
mean_h = mean(h_k);
tau_k = [tau_k{:}].'; % vettore tangente
 

%%%%%%%%%%%%%%%%%%%%%%%%%% Assembling A and F%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    y_q = h_k(j)/2 * xq_on_diag;
   
    integrand = besselh(0, 1, k*y_q);  
    A(j, j) = (1i/2)*(h_k(j)/2)*wq_on_diag.'*integrand;
end

%%%%%%%%%%%%%%%%%%%%%%% Solve the linear system %%%%%%%%%%%%%%%%%%%%%%%%%%%

psi = A\F;

%%%%%%%%%%%%%%%%%%%%%%%%%% Far field pattern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = linspace(0, 2*pi, 360)';
d = cos(theta) + 1i*sin(theta);

[xq_plot, wq_plot] = gaussquad(n_gauss_pts_farField);

farField = zeros(size(theta));

coeff = exp(1i*pi/4)/sqrt(8*pi*k);

for j = 1:N

    yq = p_k(j) + h_k(j)*xq_plot.*tau_k(j);

    integrand = exp(-1i*k*real(d*conj(yq.'))); % matrice length(d) x lenght(yq)
    farField = farField + coeff * psi(j) * h_k(j) * integrand * wq_plot;
end

minFarField = min(log10(abs(farField)));

polarplot(theta, -min(log10(abs(farField))) + log10(abs(farField)), lineWidth = 2)
disp(minFarField)

function [x, w] = gaussquad(q)
% nodi e pesi di quadratura di Gauss su [0 1]
B = ( 1:(q-1) )./ sqrt( 4*( 1:(q-1) ).^2 -1 );
[V, D] = eig( diag(B, -1) + diag(B, 1) );
x = ( diag(D)+1 )/2;
w = ( V(1, :).*V(1, :) )';
end




