function [u_scat, u_tot, times, mean_h] = BEM_func(N, k, v, u_inc)
    x_lim = [-1.5 1.5]; 
    y_lim = [-1.5 1.5];
    n_gauss_pts_plot = 5;
    n_gauss_pts_off_diag = 5;
    n_gauss_pts_on_diag = 30;
    n_points_plot = 150;

    if isa(u_inc, "double")
        u_inc=@(x) exp(1i * k * real(x*exp(-1i*u_inc)));
    end

    assert(isa(u_inc, "function_handle"))
    
    % Geometry
    v = [v v(1)];
    n_sides = length(v)-1;
    perimeter = sum(abs(v - [v(2:end) v(1)]));
    side_length = abs(diff(v));
    side_percent = side_length./perimeter;
    assert(abs( sum(side_percent)-1 ) < 1e-4);
    
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
    
    % delete the edges of Gamma (where diff(p_k)=0)
    p_k = p_k(diff(p_k) ~= 0); 
    h_k = [h_k{:}].';
    mean_h = mean(h_k);
    tau_k = [tau_k{:}].';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Assemblaggio A e F %%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % fix diag A
    [xq_on_diag, wq_on_diag] = gaussquad(n_gauss_pts_on_diag);
    for j = 1:N
        y_q = h_k(j)/2 * xq_on_diag;
       
        integrand = besselh(0, 1, k*y_q);  
        A(j, j) = (1i/2)*(h_k(j)/2)*wq_on_diag.'*integrand;
    end
    time_assembling = toc;
    
    tic
    psi = A\F;
    time_lin_sist = toc; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_plot = linspace(x_lim(1), x_lim(2), n_points_plot+1);
    y_plot = linspace(y_lim(1), y_lim(2), n_points_plot+1);
    [X, Y] = meshgrid(x_plot(1:end-1), y_plot(1:end-1)); % confronto con MPSpack
    Z = X + 1i*Y;
    
    u_scat = zeros(size(Z));
    in = inpolygon(X, Y, real(v), imag(v));
    [xq_plot, wq_plot] = gaussquad(n_gauss_pts_plot);
    xq_plot = reshape(xq_plot, 1, 1, n_gauss_pts_plot);
    wq_plot = reshape(wq_plot, 1, 1, n_gauss_pts_plot);

    for j = 1:N

        y_q = p_k(j) +  h_k(j) * xq_plot  * tau_k(j);
        r = abs(Z - y_q);

        integrand = besselh(0, 1, k*r);
        u_scat = u_scat + sum( (1i/4)  * integrand .* wq_plot * psi(j)*h_k(j), 3);
    end

    u_scat(in) = complex(NaN, NaN);
    time_plot = toc;
    u_tot = u_scat + u_inc(Z);
    u_tot(in) = complex(NaN, NaN);
    times = [time_assembling, time_lin_sist, time_plot];
end
