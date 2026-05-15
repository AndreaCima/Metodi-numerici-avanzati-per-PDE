function [psi, times, mean_h] = BEM_psi(N, k, v, theta)
    % Parameters
%     v = 0.5*[-1+1i, -1-1i, 1-1i, 1+1i]; % vertices of \Gamma (counterclockwise)
%     k = 20;
%     N = 4*k; % degrees of freedom. (Usually N ~ k ~ 1/h). Here I use N=4*k


%     theta = -pi/4;

    n_gauss_pts_off_diag = 5;
    n_gauss_pts_on_diag = 30;

    u_inc=@(x) exp(1i * k * real(x*exp(-1i*theta))); 
    
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
    tau_k = [tau_k{:}].';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Assembling A and F%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        y_q = h_k(j)/2 * xq_on_diag;
       
        integrand = besselh(0, 1, k*y_q);  
        A(j, j) = (1i/2)*(h_k(j)/2)*wq_on_diag.'*integrand;
    end
    time_assembling = toc;
    
    %%%%%%%%%%%%%%%%%%%%%%% Solve the linear system %%%%%%%%%%%%%%%%%%%%%%%
    tic
    psi = A\F;
    time_lin_sist = toc; 
    
    times = [time_assembling, time_lin_sist];
end

