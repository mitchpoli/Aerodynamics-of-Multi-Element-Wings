clc;
clear;
close all;

function V_iz_11 = calculate_self_induced_downwash_11(gamma_1,b1)
    V_iz_11 = gamma_1/(2*b1);
end

function V_iz_22 = calculate_self_induced_downwash_22(gamma_2,b2)
    V_iz_22 = gamma_2/(2*b2);
end

function epsilon_11 = calculate_self_induced_alpha_11(V_inf, gamma_01, b1)
    V_iz_11 = calculate_self_induced_downwash_11(gamma_01,b1);
    epsilon_11 = V_iz_11/V_inf;
end

function epsilon_22 = calculate_self_induced_alpha_22(V_inf, gamma_02, b2)
    V_iz_22 = calculate_self_induced_downwash_22(gamma_02,b2);
    epsilon_22 = V_iz_22/V_inf;
end

function S_t_prime = calculate_stagger(S_t, G, geometric_alpha)
    S_t_prime = S_t*cos(geometric_alpha) + G*sin(geometric_alpha);
end

function G_prime = calculate_gap(S_t, G, geometric_alpha)
    G_prime = -S_t*sin(geometric_alpha) + G*cos(geometric_alpha);
end

function r_distance = calculate_distance(S_t, G, geometric_alpha, delta_y)
    S_t_prime = calculate_stagger(S_t, G, geometric_alpha);
    G_prime =  calculate_gap(S_t, G, geometric_alpha);
    r_distance = sqrt(delta_y^2+S_t_prime^2+G_prime^2);
end

function sin_phi = calculate_sin_phi(S_t, G, geometric_alpha, delta_y)
    S_t_prime = calculate_stagger(S_t, G, geometric_alpha);
    G_prime =  calculate_gap(S_t, G, geometric_alpha);
    r_distance = calculate_distance(S_t, G, geometric_alpha, delta_y);
    sin_phi = sqrt(S_t_prime^2+G_prime^2)/r_distance;
end

function cos_phi = calculate_cos_phi(S_t, G, geometric_alpha)
    S_t_prime = calculate_stagger(S_t, G, geometric_alpha);
    G_prime =  calculate_gap(S_t, G, geometric_alpha);
    cos_phi = S_t_prime/sqrt(S_t_prime^2+G_prime^2);
end

function r_prime_val = calculate_r_prime(S_t, G, geometric_alpha, delta_y)
    % Calculate perpendicular distance to trailing vortex filament
    % r' is the perpendicular distance from point to trailing vortex filament
    G_prime = calculate_gap(S_t, G, geometric_alpha);
    r_prime_val = sqrt(G_prime^2 + delta_y^2);
end

function V_iz_12_1 = calculate_bound_vortex_speed_12_1(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_02)

    V_iz_12_1 = zeros(1, N_control_points-1);
    
    % θ per i punti di controllo dell'ala 1 (come prima)
    theta_1_range = linspace(0, pi, N_control_points+1);
    % θ di integrazione lungo l'ala 2
    theta_2_edges = linspace(0, pi, N_integral_points+1);
    theta_2_range = 0.5*(theta_2_edges(1:end-1) + theta_2_edges(2:end));

    
    cos_phi = calculate_cos_phi(S_t, G, geometric_alpha);

    for i = 1:N_control_points-1
        % Punto di controllo i sull'ala 1
        theta_1 = (theta_1_range(i) + theta_1_range(i+1)) / 2;
        y_1 = -b1/2 * cos(theta_1);

        integral_sum = 0;

        for j = 1:N_integral_points
            theta_2 = theta_2_range(j);

            % Posizione della sezione j sull'ala 2
            y_2 = -b2/2 * cos(theta_2);
            delta_y = y_2 - y_1;

            % Distanza e geometria
            r = calculate_distance(S_t, G, geometric_alpha, delta_y);
            sin_phi = calculate_sin_phi(S_t, G, geometric_alpha, delta_y);

            % *** Γ₂(θ₂) interpolata dalla distribuzione gamma_02 ***
            Gamma_2 = interpolate_gamma(theta_2, gamma_02);

            % Jacobiano dy/dθ = (b2/2) sin θ₂
            dy_dtheta = (b2/2) * sin(theta_2);

            % Kernel di Biot–Savart per il bound vortex
            integrand = (Gamma_2 * sin_phi * dy_dtheta) / (4*pi*r^2);

            dtheta = pi / N_integral_points;
            integral_sum = integral_sum + integrand * dtheta;
        end

        % Proiezione nella direzione z
        V_iz_12_1(i) = -cos_phi * integral_sum;
    end
end

function V_iz_21_1 = calculate_bound_vortex_speed_21_1(N_control_points, N_integral_points, b2, geometric_alpha, S_t, G, b1, gamma_01)

    V_iz_21_1 = zeros(1, N_control_points-1);
    
    theta_2_range = linspace(0, pi, N_control_points+1);
    theta_1_range = linspace(0, pi, N_integral_points+1);
    
    cos_phi = calculate_cos_phi(S_t, G, geometric_alpha);

    for i = 1:N_control_points-1
        theta_2 = (theta_2_range(i) + theta_2_range(i+1)) / 2;
        y_2 = -b2/2 * cos(theta_2);

        integral_sum = 0;

        for j = 1:N_integral_points
            theta_1 = theta_1_range(j);
            y_1 = -b1/2 * cos(theta_1);
            delta_y = y_2 - y_1;

            r = calculate_distance(S_t, G, geometric_alpha, delta_y);
            sin_phi = calculate_sin_phi(S_t, G, geometric_alpha, delta_y);

            % *** Γ₁(θ₁) interpolata ***
            Gamma_1 = interpolate_gamma(theta_1, gamma_01);

            dy_dtheta = (b1/2) * sin(theta_1);

            integrand = (Gamma_1 * sin_phi * dy_dtheta) / (4*pi*r^2);

            dtheta = pi / N_integral_points;
            integral_sum = integral_sum + integrand * dtheta;
        end

        V_iz_21_1(i) = cos_phi * integral_sum;
    end
end

function V_iz_12_2 = calculate_trailing_vortex_speed_12_2(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_02)

    V_iz_12_2 = zeros(1, N_control_points-1);
    theta_1_range = linspace(0, pi, N_control_points+1);
    theta_2_range = linspace(0, pi, N_integral_points+1);
    S_t_prime = calculate_stagger(S_t, G, geometric_alpha);

    for i = 1:N_control_points-1
        theta_1 = (theta_1_range(i)+theta_1_range(i+1))/2;
        y_1 = -b1/2 * cos(theta_1);

        integral_sum = 0;

        for j = 1:N_integral_points
            theta_2 = theta_2_range(j);
            y_2 = -b2/2 * cos(theta_2);

            delta_y_2_1 = y_2 - y_1;

            r_prime = calculate_r_prime(S_t, G, geometric_alpha, delta_y_2_1);
            first_param = delta_y_2_1 / r_prime^2;

            r = calculate_distance(S_t, G, geometric_alpha, delta_y_2_1);
            second_param = 1 - S_t_prime / r;

            % *** Γ₂(θ₂) interpolata, senza cos(θ2) ***
            Gamma_2 = interpolate_gamma(theta_2, gamma_02);

            integrand = first_param * second_param * Gamma_2;

            dtheta = pi / N_integral_points;
            integral_sum = integral_sum + integrand * dtheta;
        end

        V_iz_12_2(i) = integral_sum / (4*pi);
    end
end

function V_iz_21_2 = calculate_trailing_vortex_speed_21_2(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_01)
    V_iz_21_2 = zeros(1, N_control_points-1);
    theta_2_range = linspace(0, pi, N_control_points+1);
    theta_1_range = linspace(0, pi, N_integral_points+1);
    S_t_prime = calculate_stagger(S_t, G, geometric_alpha);

    for i = 1:N_control_points-1
        theta_2 = (theta_2_range(i)+theta_2_range(i+1))/2;
        y_2 = -b2/2 * cos(theta_2);

        integral_sum = 0;

        for j = 1:N_integral_points
            theta_1 = theta_1_range(j);
            y_1 = -b1/2 * cos(theta_1);
            delta_y_2_1 = y_2-y_1;
            r_prime = calculate_r_prime(S_t, G, geometric_alpha, delta_y_2_1);
            first_param = delta_y_2_1/r_prime^2;
           
            r = calculate_distance(S_t, G, geometric_alpha, delta_y_2_1);
            second_param = 1 + S_t_prime/r;
            
            % *** Γ1(θ1) interpolata, senza cos(θ2) ***
            Gamma_1 = interpolate_gamma(theta_1, gamma_01);

            integrand = first_param*second_param*Gamma_1;

            dtheta = pi/N_integral_points;

            integral_sum = integral_sum + integrand * dtheta;

        end
        V_iz_21_2(i) = integral_sum/(4*pi);
    end
end

function V_iz_12 = calculate_V_iz_12(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_02)
    V_iz_12_1 = calculate_bound_vortex_speed_12_1(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_02);   
    V_iz_12_2 = calculate_trailing_vortex_speed_12_2(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_02);
    V_iz_12 = V_iz_12_2+V_iz_12_1;
end

function V_iz_21 = calculate_V_iz_21(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_01)
    V_iz_21_1 = calculate_bound_vortex_speed_21_1(N_control_points, N_integral_points, b2, geometric_alpha,S_t,G,b1,gamma_01);   
    V_iz_21_2 = calculate_trailing_vortex_speed_21_2(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_01);
    V_iz_21 = V_iz_21_1+V_iz_21_2;
end

function epsilon_1_2 = calculate_induced_angle_1_2(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_02,V_inf)
    V_iz_12 = calculate_V_iz_12(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_02);
    epsilon_1_2 = V_iz_12/V_inf;
end

function epsilon_2_1 = calculate_induced_angle_2_1(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_01,V_inf)
    V_iz_21 =calculate_V_iz_21(N_control_points, N_integral_points, b1, geometric_alpha,S_t,G,b2,gamma_01);
    epsilon_2_1 = V_iz_21/V_inf;
end

function alpha_eff_1 = calculate_alpha_eff_1(geometric_alpha, N_control_points, N_integral_points, b1, S_t, G, b2, gamma_02, V_inf, gamma_01)

    % Self-induced: scalare
    epsilon_11 = calculate_self_induced_alpha_11(V_inf, gamma_01, b1);

    % Mutual-induced: vettore
    epsilon_12 = calculate_induced_angle_1_2(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_02, V_inf);

    alpha_eff_1 = geometric_alpha - epsilon_11 - epsilon_12;
end

function alpha_eff_2 = calculate_alpha_eff_2(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_01,V_inf, gamma_02)
    epsilon_22 = calculate_self_induced_alpha_22(V_inf, gamma_02, b2);
    epsilon_21 = calculate_induced_angle_2_1(N_control_points, N_integral_points, b2, geometric_alpha, S_t, G, b1, gamma_01, V_inf);

    alpha_eff_2 = geometric_alpha - epsilon_22 - epsilon_21;

end

function [gamma_01_initial, gamma_02_initial] = initial_gamma_guess(V_inf, b1, b2, S1, S2, N_control_points, alpha, model)

    % Generate theta at panel MIDPOINTS (length = N_control_points - 1)
    theta_edges = linspace(0, pi, N_control_points);
    theta = 0.5*(theta_edges(1:end-1) + theta_edges(2:end));

    switch model
        
        case 'elliptic'
            % OLD elliptical guess
            C_L_isolated = 2*pi*alpha;
            A1_1 = (2 * V_inf * S1 * C_L_isolated)/(pi*b1);
            A1_2 = (2 * V_inf * S2 * C_L_isolated)/(pi*b2);

            gamma_01_initial = A1_1 * sin(theta);
            gamma_02_initial = A1_2 * sin(theta);

        case 'gaussian'
            % Gaussian guess
            x = linspace(-1, 1, N_control_points - 1);

            gamma_02_initial = -0.010 * exp(-460.52 * x.^2) + 0.217;
            gamma_01_initial =  0.003 * exp(-460.52 * x.^2) + 0.237;

        case 'sinusoidal'
            % sum of sinusoids
            C_L_isolated = 2*pi*alpha;

            A1_1 = (2 * V_inf * S1 * C_L_isolated)/(pi*b1); 
            A1_2 = (2 * V_inf * S2 * C_L_isolated)/(pi*b2);

            % Introduce harmonics for non-elliptic shapes
            A2_1 = 0.20 * A1_1;
            A3_1 = -0.10 * A1_1;

            A2_2 = -0.15 * A1_2;
            A3_2 =  0.10 * A1_2;

            gamma_01_initial = A1_1*sin(theta) + A2_1*sin(2*theta) + A3_1*sin(3*theta);
            gamma_02_initial = A1_2*sin(theta) + A2_2*sin(2*theta) + A3_2*sin(3*theta);

    end

end

function gamma_1 = calculate_gamma_1(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_02,V_inf, gamma_01,c)
    alpha_eff_1 = calculate_alpha_eff_1(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_02,V_inf, gamma_01);
    gamma_1 = zeros(1, length(alpha_eff_1));  % Match epsilon_12 size
    for i = 1:length(alpha_eff_1)
        gamma_1(i) = 0.5 * V_inf * c * 2 * pi * alpha_eff_1(i);
    end
end

function gamma_2 = calculate_gamma_2(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_02,V_inf, gamma_01,c)
    alpha_eff_2 = calculate_alpha_eff_2(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_01,V_inf, gamma_02);
    gamma_2 = zeros(1, length(alpha_eff_2));  % Match epsilon_12 size
    for i = 1:length(alpha_eff_2)
       gamma_2(i) = 0.5 * V_inf * c * 2 * pi * alpha_eff_2(i);
    end
end

function [lift_coefficient_1, lift_coefficient_2] = solve_my_problem()
    % Main solver for tandem wing lifting-line problem
    
    % ========== INPUT PARAMETERS ==========
    % Flow conditions
    V_inf = 10;                    % Freestream velocity [m/s]
    geometric_alpha = 5 * pi/180;     % Geometric angle of attack [rad]
    rho = 1.2;
    
    % Wing 1 (Fore wing) parameters
    b1 = 1.0;                         % Span [m]
    S1 = 0.1;                         % Wing area [m²]
    
    % Wing 2 (Hind wing) parameters  
    b2 = 1.0;                         % Span [m]
    S2 = 0.1;                         % Wing area [m²]

    % Chord
    c = 0.1;                         % Chord [m]
    
    % Tandem configuration
    S_t = 0;                  % Stagger [m]
    G = 10 * c;                    % Gap [m] (negative = hind above fore)
    
    % Numerical parameters
    N_control_points = 100;            % Control points per wing
    N_integral_points = 1000;          % Integration points
    max_iterations = 2;             % Maximum iterations
    tolerance = 1e-3;                 % Convergence tolerance

    
    [gamma_01_initial, gamma_02_initial] = initial_gamma_guess(V_inf, b1, b2, S1, S2, N_control_points, geometric_alpha, 'sinusoidal');
    gamma_1_prev = gamma_01_initial;
    gamma_2_prev = gamma_02_initial;
    
    for iter = 1:max_iterations
        gamma_1_new = calculate_gamma_1(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_2_prev,V_inf, gamma_1_prev,c);
        gamma_2_new = calculate_gamma_2(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_2_prev,V_inf, gamma_1_prev,c);

       % Check convergence
       diff_1 = max(abs(gamma_1_new - gamma_1_prev));
       diff_2 = max(abs(gamma_2_new - gamma_2_prev));
       max_diff = max(diff_1, diff_2);
        
       fprintf('  Max change: %.6f\n', max_diff);
        
       if max_diff < tolerance
           fprintf('Converged after %d iterations!\n', iter);
           break;
       end

       % Update for next iteration (with relaxation if needed)
       relaxation = 0.1;  % Under-relaxation factor
       gamma_1_prev = relaxation * gamma_1_new + (1 - relaxation) * gamma_1_prev;
       gamma_2_prev = relaxation * gamma_2_new + (1 - relaxation) * gamma_2_prev;

       if iter == max_iterations
        fprintf('Reached maximum iterations (%d)\n', max_iterations);
       end
    end
    % ========== OUTPUT RESULTS ==========
    gamma_1_final = gamma_1_prev;
    gamma_2_final = gamma_2_prev;

    % Calculate total circulation (integral along span)
    total_gamma_1_final = sum(gamma_1_final) * (b1/(N_control_points-1));  % Approximate integral
    total_gamma_2_final = sum(gamma_2_final) * (b2/(N_control_points-1));  % FIXED: gamma_2_final

    lift_1 = rho * V_inf * total_gamma_1_final;
    lift_2 = rho * V_inf * total_gamma_2_final;

    lift_coefficient_1 = lift_1 / (0.5 * rho * V_inf^2 * S1);  
    lift_coefficient_2 = lift_2 / (0.5 * rho * V_inf^2 * S2);  
end

function [lift_coefficient_1, lift_coefficient_2, debug_info] = solve_my_problem_debug()
    % Debug solver for tandem wing lifting-line problem
    
    % ========== INPUT PARAMETERS ==========
    % Flow conditions
    V_inf = 14;                    % Freestream velocity [m/s]
    geometric_alpha = 5 * pi/180;     % Geometric angle of attack [rad]
    rho = 1.225;
    
    % Wing 1 (Fore wing) parameters
    b1 = 1.0;                         % Span [m]
    S1 = 0.1;                         % Wing area [m²]
    
    % Wing 2 (Hind wing) parameters  
    b2 = 1.0;                         % Span [m]
    S2 = 0.1;                         % Wing area [m²]

    % Chord
    c = 0.1;                         % Chord [m]
    
    % Tandem configuration
    S_t = 1.63 * c;                  % Stagger [m]
    G = -0.5 * c;                    % Gap [m] (negative = hind above fore)
    
    % Numerical parameters
    N_control_points = 40;            % Control points per wing
    N_integral_points = 100;          % Integration points
    max_iterations = 100;             % Maximum iterations
    tolerance = 1e-6;                 % Convergence tolerance

    
    [gamma_01_initial, gamma_02_initial] = initial_gamma_guess(V_inf, b1, b2, S1, S2, N_control_points, geometric_alpha, 'sinusoidal');
    gamma_1_prev = gamma_01_initial;
    gamma_2_prev = gamma_02_initial;
    
    fprintf('=== STARTING SOLVER DEBUG ===\n');
    
    for iter = 1:max_iterations
        gamma_1_new = calculate_gamma_1(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_2_prev,V_inf, gamma_1_prev,c);
        gamma_2_new = calculate_gamma_2(geometric_alpha, N_control_points, N_integral_points, b1,S_t,G,b2,gamma_2_prev,V_inf, gamma_1_prev,c);

       % Check convergence
       diff_1 = max(abs(gamma_1_new - gamma_1_prev));
       diff_2 = max(abs(gamma_2_new - gamma_2_prev));
       max_diff = max(diff_1, diff_2);
        
       fprintf('Iter %d: Max change: %.6f\n', iter, max_diff);
        
       if max_diff < tolerance
           fprintf('=== CONVERGED AFTER %d ITERATIONS ===\n', iter);
           
           % Calculate and print all debug information
           print_debug_info(geometric_alpha, N_control_points, N_integral_points, ...
                           b1, S_t, G, b2, V_inf, c, gamma_1_new, gamma_2_new);
           break;
       end

       % Update for next iteration
       relaxation = 0.2;  % Under-relaxation factor
       gamma_1_prev = relaxation * gamma_1_new + (1 - relaxation) * gamma_1_prev;
       gamma_2_prev = relaxation * gamma_2_new + (1 - relaxation) * gamma_2_prev;

       if iter == max_iterations
        fprintf('Reached maximum iterations (%d)\n', max_iterations);
        % Print debug info even if not fully converged
        print_debug_info(geometric_alpha, N_control_points, N_integral_points, ...
                        b1, S_t, G, b2, V_inf, c, gamma_1_prev, gamma_2_prev);
       end
    end
    
    % ========== OUTPUT RESULTS ==========
    gamma_1_final = gamma_1_prev;
    gamma_2_final = gamma_2_prev;

    % Calculate total circulation (integral along span)
    total_gamma_1_final = sum(gamma_1_final) * (b1/(N_control_points-1));
    total_gamma_2_final = sum(gamma_2_final) * (b2/(N_control_points-1));

    lift_1 = rho * V_inf * total_gamma_1_final;
    lift_2 = rho * V_inf * total_gamma_2_final;

    lift_coefficient_1 = lift_1 / (0.5 * rho * V_inf^2 * S1);  
    lift_coefficient_2 = lift_2 / (0.5 * rho * V_inf^2 * S2);
    
    % Store debug info
    debug_info.CL1 = lift_coefficient_1;
    debug_info.CL2 = lift_coefficient_2;
    debug_info.gamma_1 = gamma_1_final;
    debug_info.gamma_2 = gamma_2_final;
    debug_info.iterations = iter;
end

function print_debug_info(geometric_alpha, N_control_points, N_integral_points, b1, S_t, G, b2, V_inf, c, gamma_1, gamma_2)
    fprintf('\n=== DEBUG INFORMATION ===\n');
    
    % Calculate all intermediate values
    epsilon_11 = calculate_self_induced_alpha_11(V_inf, gamma_1(1), b1);  % Use first element for scalar
    epsilon_22 = calculate_self_induced_alpha_22(V_inf, gamma_2(1), b2);
    
    V_iz_12 = calculate_V_iz_12(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_2);
    V_iz_21 = calculate_V_iz_21(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_1);
    
    epsilon_12 = V_iz_12 / V_inf;
    epsilon_21 = V_iz_21 / V_inf;
    
    alpha_eff_1 = calculate_alpha_eff_1(geometric_alpha, N_control_points, N_integral_points, b1, S_t, G, b2, gamma_2, V_inf, gamma_1);
    alpha_eff_2 = calculate_alpha_eff_2(geometric_alpha, N_control_points, N_integral_points, b1, S_t, G, b2, gamma_1, V_inf, gamma_2);
    
    % Print summary at root section (first control point)
    fprintf('\n--- ROOT SECTION (Mid-span) VALUES ---\n');
    fprintf('Geometric alpha: %.4f°\n', geometric_alpha * 180/pi);
    fprintf('Self-induced angles:\n');
    fprintf('  epsilon_11 (fore wing): %.6f°\n', epsilon_11 * 180/pi);
    fprintf('  epsilon_22 (hind wing): %.6f°\n', epsilon_22 * 180/pi);
    fprintf('\nMutual induced velocities [m/s]:\n');
    fprintf('  V_iz_12 (hind→fore): %.6f\n', V_iz_12(1));
    fprintf('  V_iz_21 (fore→hind): %.6f\n', V_iz_21(1));
    fprintf('\nMutual induced angles:\n');
    fprintf('  epsilon_12 (hind→fore): %.6f°\n', epsilon_12(1) * 180/pi);
    fprintf('  epsilon_21 (fore→hind): %.6f°\n', epsilon_21(1) * 180/pi);
    fprintf('\nEffective angles of attack:\n');
    fprintf('  alpha_eff_1 (fore): %.4f°\n', alpha_eff_1(1) * 180/pi);
    fprintf('  alpha_eff_2 (hind): %.4f°\n', alpha_eff_2(1) * 180/pi);
    fprintf('\nCirculations [m²/s]:\n');
    fprintf('  gamma_1 (fore): %.6f\n', gamma_1(1));
    fprintf('  gamma_2 (hind): %.6f\n', gamma_2(1));
    
    % Print spanwise distribution (first 5 points)
    fprintf('\n--- SPANWISE DISTRIBUTION (First 5 points) ---\n');
    fprintf('Station  Gamma_1    Gamma_2    alpha_eff_1  alpha_eff_2\n');
    for i = 1:min(5, length(gamma_1))
        fprintf('%4d    %8.4f   %8.4f   %8.4f°    %8.4f°\n', ...
                i, gamma_1(i), gamma_2(i), alpha_eff_1(i)*180/pi, alpha_eff_2(i)*180/pi);
    end
    
    % Print min/max values
    fprintf('\n--- EXTREMUM VALUES ---\n');
    fprintf('Gamma_1: min=%.4f, max=%.4f, mean=%.4f\n', min(gamma_1), max(gamma_1), mean(gamma_1));
    fprintf('Gamma_2: min=%.4f, max=%.4f, mean=%.4f\n', min(gamma_2), max(gamma_2), mean(gamma_2));
    fprintf('Alpha_eff_1: min=%.2f°, max=%.2f°\n', min(alpha_eff_1)*180/pi, max(alpha_eff_1)*180/pi);
    fprintf('Alpha_eff_2: min=%.2f°, max=%.2f°\n', min(alpha_eff_2)*180/pi, max(alpha_eff_2)*180/pi);
    
    % Print induced velocity components
    V_iz_12_1 = calculate_bound_vortex_speed_12_1(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_2);
    V_iz_12_2 = calculate_trailing_vortex_speed_12_2(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_2);
    V_iz_21_1 = calculate_bound_vortex_speed_21_1(N_control_points, N_integral_points, b2, geometric_alpha, S_t, G, b1, gamma_1);
    V_iz_21_2 = calculate_trailing_vortex_speed_21_2(N_control_points, N_integral_points, b1, geometric_alpha, S_t, G, b2, gamma_1);
    
    fprintf('\n--- INDUCED VELOCITY COMPONENTS [m/s] ---\n');
    fprintf('V_iz_12_1 (bound): %.6f\n', V_iz_12_1(1));
    fprintf('V_iz_12_2 (trailing): %.6f\n', V_iz_12_2(1));
    fprintf('V_iz_21_1 (bound): %.6f\n', V_iz_21_1(1));
    fprintf('V_iz_21_2 (trailing): %.6f\n', V_iz_21_2(1));
    
    fprintf('=== END DEBUG INFORMATION ===\n\n');
end

function [gaps, total_CLs, CL1s, CL2s, lift_ratios] = solve_gap_sweep(gap_range)
    % Solver for tandem wing that sweeps over gap values
    % Input: gap_range - array of gap values to test [m]
    % Output: gaps, total_CLs, CL1s, CL2s, lift_ratios - results for plotting
    
    % ========== FIXED PARAMETERS ==========
    % Flow conditions
    V_inf = 100;                       % Freestream velocity [m/s]
    geometric_alpha = 10 * pi/180;     % Geometric angle of attack [rad]
    rho = 1.2;
    
    % Wing 1 (Fore wing) parameters
    b1 = 0.4;                         % Span [m]
    S1 = 0.04;                         % Wing area [m²]
    
    % Wing 2 (Hind wing) parameters  
    b2 = 0.4;                         % Span [m]
    S2 = 0.04;                         % Wing area [m²]

    % Chord
    c = 0.1;                         % Chord [m]
    
    % Tandem configuration (gap varies, stagger fixed)
    S_t = 0.5 * c;                  % Stagger [m] (constant)
    
    % Numerical parameters
    N_control_points = 40;            % Control points per wing
    N_integral_points = 200;          % Integration points
    max_iterations = 20;              % Maximum iterations
    tolerance = 1e-3;                 % Convergence tolerance

    % Initialize results arrays
    gaps = gap_range;
    total_CLs = zeros(size(gap_range));
    CL1s = zeros(size(gap_range));
    CL2s = zeros(size(gap_range));
    lift_ratios = zeros(size(gap_range));
    
    fprintf('=== GAP SWEEP ANALYSIS ===\n');
    fprintf('Testing %d gap values from %.2fc to %.2fc\n', ...
            length(gap_range), min(gap_range)/c, max(gap_range)/c);
    
    % Calculate single wing lift coefficient for comparison
    single_wing_CL = 2 * pi * geometric_alpha;
    fprintf('Single wing CL (2π×α): %.4f\n', single_wing_CL);
    
    % Loop over all gap values
    for k = 1:length(gap_range)
        G = -gaps(k);
        
        fprintf('\n--- Gap = %.2fc (%.3f m) ---\n', G/c, G);
        
        % Initialize for this gap value
        [gamma_01_initial, gamma_02_initial] = initial_gamma_guess(V_inf, b1, b2, S1, S2, N_control_points, geometric_alpha, 'sinusoidal');
        gamma_1_prev = gamma_01_initial;
        gamma_2_prev = gamma_02_initial;
        
        % Iterative solver for this gap
        for iter = 1:max_iterations
            gamma_1_new = calculate_gamma_1(geometric_alpha, N_control_points, N_integral_points, b1, S_t, G, b2, gamma_2_prev, V_inf, gamma_1_prev, c);
            gamma_2_new = calculate_gamma_2(geometric_alpha, N_control_points, N_integral_points, b1, S_t, G, b2, gamma_2_prev, V_inf, gamma_1_prev, c);

            % Check convergence
            diff_1 = max(abs(gamma_1_new - gamma_1_prev));
            diff_2 = max(abs(gamma_2_new - gamma_2_prev));
            max_diff = max(diff_1, diff_2);
            
            if max_diff < tolerance
                fprintf('  Converged after %d iterations\n', iter);
                break;
            end

            % Update for next iteration
            relaxation = 0.3;
            gamma_1_prev = relaxation * gamma_1_new + (1 - relaxation) * gamma_1_prev;
            gamma_2_prev = relaxation * gamma_2_new + (1 - relaxation) * gamma_2_prev;

            if iter == max_iterations
                fprintf('  Reached max iterations (%d), diff: %.6f\n', max_iterations, max_diff);
            end
        end
        
        % Calculate results for this gap
        gamma_1_final = gamma_1_prev;
        gamma_2_final = gamma_2_prev;

        % Calculate total circulation
        theta_1_range = linspace(0, pi, N_control_points+1);  % N+1 elementi (bordi pannelli)
        % I punti di controllo sono ai CENTRI dei pannelli
        theta_1_ctrl = zeros(1, N_control_points-1);
        for i = 1:N_control_points-1
            theta_1_ctrl(i) = (theta_1_range(i) + theta_1_range(i+1)) / 2;
        end
        y_1 = -(b1/2) * cos(theta_1_ctrl);  % Ora: length(y_1) = N_control_points-1

        % Integrazione corretta
        total_gamma_integral_1 = trapz(y_1, gamma_1_final);  % Entrambi hanno length N_control_points-1

        CL1 = (2 * total_gamma_integral_1) / (V_inf * S1);

        % Stesso per Wing 2:
        theta_2_range = linspace(0, pi, N_control_points+1);
        theta_2_ctrl = zeros(1, N_control_points-1);
        for i = 1:N_control_points-1
            theta_2_ctrl(i) = (theta_2_range(i) + theta_2_range(i+1)) / 2;
        end
        y_2 = -(b2/2) * cos(theta_2_ctrl);

        total_gamma_integral_2 = trapz(y_2, gamma_2_final);

        CL2 = (2 * total_gamma_integral_2) / (V_inf * S2);

        total_CL = (CL1 * S1 + CL2 * S2) / (S1 + S2);
        
        % Calculate lift ratio (tandem efficiency vs single wing)
        lift_ratio = total_CL / single_wing_CL;
        
        % Store results
        CL1s(k) = CL1;
        CL2s(k) = CL2;
        total_CLs(k) = total_CL;
        lift_ratios(k) = lift_ratio;
        
        fprintf('  Results: CL1=%.4f, CL2=%.4f, Total_CL=%.4f, Lift_Ratio=%.4f\n', CL1, CL2, total_CL, lift_ratio);
    end
    
    % Plot results
    plot_gap_sweep_results(gaps/c, CL1s, CL2s, total_CLs, lift_ratios, single_wing_CL, geometric_alpha);
    
    fprintf('\n=== GAP SWEEP COMPLETE ===\n');
end
%% 
function plot_gap_sweep_results(gaps_chord, CL1s, CL2s, total_CLs, lift_ratios, single_wing_CL, geometric_alpha)
    % Plot the gap sweep results with connected lines
    
    % Figure 1: Lift coefficients vs gap
    figure;
    
    % Plot individual and total lift coefficients with connected lines
    plot(gaps_chord, CL1s, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Fore Wing CL');
    hold on;
    plot(gaps_chord, CL2s, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Hind Wing CL');
    plot(gaps_chord, total_CLs, 'k-^', 'LineWidth', 3, 'MarkerSize', 8, 'DisplayName', 'Total CL');
    
    % Add reference line for single wing performance
    plot([min(gaps_chord), max(gaps_chord)], [single_wing_CL, single_wing_CL], 'g--', 'LineWidth', 2, ...
         'DisplayName', sprintf('Single Wing CL (%.3f)', single_wing_CL));
    
    % Formatting
    xlabel('Gap (c)');
    ylabel('Lift Coefficient C_L');
    title('Tandem Wing Lift vs Gap');
    legend('Location', 'best');
    grid on;
    
    % Find optimal gap
    [max_CL, max_idx] = max(total_CLs);
    optimal_gap = gaps_chord(max_idx);
    
    % Add optimal point annotation
    hold on;
    plot(optimal_gap, max_CL, 'mo', 'MarkerSize', 12, 'LineWidth', 3, ...
         'DisplayName', sprintf('Optimal: %.2fc', optimal_gap));
    
    fprintf('\nOptimal gap: %.2fc with Total CL = %.4f\n', optimal_gap, max_CL);
    
    % Add labels for positive/negative gap regions
    x_limits = xlim;
    y_limits = ylim;
    text(x_limits(1) + 0.1*diff(x_limits), y_limits(2) - 0.1*diff(y_limits), ...
         'Negative gap = hind wing above fore wing', 'FontSize', 10, 'Color', 'blue');
    text(x_limits(1) + 0.1*diff(x_limits), y_limits(2) - 0.15*diff(y_limits), ...
         'Positive gap = hind wing below fore wing', 'FontSize', 10, 'Color', 'red');
    
    % Figure 2: Lift ratio vs gap
    figure;
    
    % Plot lift ratio with connected line
    plot(gaps_chord, lift_ratios, 'b-o', 'LineWidth', 3, 'MarkerSize', 8);
    hold on;
    
    % Add reference line at 1.0 (equal to single wing performance)
    plot([min(gaps_chord), max(gaps_chord)], [1.0, 1.0], 'r--', 'LineWidth', 2, ...
         'DisplayName', 'Single Wing Reference');
    
    % Find optimal lift ratio
    [max_ratio, max_ratio_idx] = max(lift_ratios);
    optimal_ratio_gap = gaps_chord(max_ratio_idx);
    
    % Add optimal point annotation
    plot(optimal_ratio_gap, max_ratio, 'go', 'MarkerSize', 12, 'LineWidth', 3, ...
         'DisplayName', sprintf('Optimal: %.2fc', optimal_ratio_gap));
    
    % Formatting
    xlabel('Gap (c)');
    ylabel('Lift Ratio (Total CL / Single Wing CL)');
    title('Tandem Wing Efficiency vs Gap');
    legend('Lift Ratio', 'Single Wing Reference', 'Optimal Gap', 'Location', 'best');
    grid on;
    
    % Add interpretation labels
    x_limits = xlim;
    y_limits = ylim;
    text(x_limits(1) + 0.1*diff(x_limits), y_limits(1) + 0.9*diff(y_limits), ...
         sprintf('Max efficiency: %.1f%% at gap = %.2fc', max_ratio*100, optimal_ratio_gap), ...
         'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    
    text(x_limits(1) + 0.1*diff(x_limits), y_limits(1) + 0.8*diff(y_limits), ...
         'Ratio > 1: Tandem outperforms single wing', ...
         'FontSize', 10, 'Color', 'green', 'BackgroundColor', 'white');
    
    text(x_limits(1) + 0.1*diff(x_limits), y_limits(1) + 0.7*diff(y_limits), ...
         'Ratio < 1: Single wing outperforms tandem', ...
         'FontSize', 10, 'Color', 'red', 'BackgroundColor', 'white');

     text(x_limits(1) + 0.1*diff(x_limits), y_limits(1) + 0.6*diff(y_limits), ...
         sprintf('Angle = %.0f', geometric_alpha*180/pi), ...
         'FontSize', 10, 'Color', 'red', 'BackgroundColor', 'white');
    
    fprintf('Maximum lift ratio: %.4f at gap = %.2fc (%.1f%% efficiency)\n', ...
            max_ratio, optimal_ratio_gap, max_ratio*100);
end

%%
% Define gap range (in meters or chord multiples)
gap_range = linspace(-2*0.1, 2*0.1, 50);  % From -2c to 2c in 20 steps

% Run the gap sweep
[gaps, total_CLs, CL1s, CL2s, lift_ratios] = solve_gap_sweep(gap_range);














function Gamma_val = interpolate_gamma(theta, gamma_vec)
    % gamma_vec: valori di Γ ai punti di controllo
    % theta: angolo in cui vuoi Γ(θ)

    % Numero di punti di controllo effettivi
    N_ctrl = numel(gamma_vec);

    % Costruisci le θ alle estremità "di pannello"
    theta_range = linspace(0, pi, N_ctrl + 1);

    % Punti di controllo = punti medi dei pannelli
    theta_mid = 0.5 * (theta_range(1:end-1) + theta_range(2:end));
    % -> theta_mid ha lunghezza N_ctrl, come gamma_vec

    % Interpolazione lineare di Γ(θ)
    Gamma_val = interp1(theta_mid, gamma_vec, theta, 'linear', 'extrap');
end