 clc;
 clear; 
 close all; 

%% Import data 
load('xy.mat');
load('low10.mat');
load('low15.mat');
load('low20.mat');
load('low25.mat');
load('low30.mat');
load('up10.mat');
load('up15.mat');
load('up20.mat');
load('up25.mat');
load('up30.mat');


%% INPUT EXPERIMENTAL DATA
% xy are the values of the experimental data extrapolatex from Figure 3 (CL vs alpha) using the provided file image_digitizer.m
alpha_exp_deg = xy(:,1);
CL_exp = xy(:,2);


%% Parameters
AR = 4;     % Aspect ratio [-]
c0 = 0.1;   % Chord [m]
b = 0.4;    % Span [m]


%% Linear fit to obtain a_exp and alphaL0 from the paper results
% The results extracted in this section will be used only for comparison
% with our proper results

% Select only the almost linear part
mask = (alpha_exp_deg >= -10) & (alpha_exp_deg <= 10);

alpha_fit_deg = alpha_exp_deg(mask);
CL_fit = CL_exp(mask);

% Degrees converted to radians
alpha_fit_rad = deg2rad(alpha_fit_deg);

% Linear fit
p = polyfit(alpha_fit_rad, CL_fit, 1);

a_exp = p(1);                 % Experimental slope [1/rad]
alphaL0 = -p(2)/p(1);         % Angle at CL = 0 [rad]
alphaL0_deg = rad2deg(alphaL0);

fprintf('=== Experimental fit CL(alpha) ===\n');
fprintf('Slope from experimental data = %.3f [1/rad]\n', a_exp);
fprintf('alphaL0 = %.2f [deg]\n', alphaL0_deg);



%% Profile parameter for thin airfoil theory
% Values are extrapolatex from Figure 4 using the provided file image_digitizer.m
% In this section, the coordinates of the profile where there is no
% detachment are extrapolated. In addition, the maximum chamber M and its position P are calculated.

smooth_win = 6;
Npts = 1000;

P_real = 0;
M_real = 0;

% xu: x coordinate of the upper profile
% yu: y coordinate of the upper profile
% xl: x coordinate of the lower profile
% yl: y coordinate of the lower profile
% upi: coordinate extrapolated from the profile at i degree of the upper profile
% lowi: coordinate extrapolated from the profile at i degree of the lower profile

% Cell arrays to store camber lines and grids
xc_cells = cell(1,6);
yc_cells = cell(1,6);
c_values = zeros(1,6);   % chord values for each extracted profile

% Case 1: flat plate
xc_cells{1} = linspace(0,1,Npts)';
yc_cells{1} = zeros(Npts,1);
c_values(1) = c0;

% List of imported datasets
up_list  = {[], up10, up15, up20, up25, up30};
low_list = {[], low10, low15, low20, low25, low30};

% Preallocate M and P arrays
M = zeros(1,6);
P = zeros(1,6);

% Loop over extracted profiles (10°, 15°, 20°, 25°, 30°)
for i = 2:6
    
    % Upper and lower surface datasets for case i
    up_data  = up_list{i};
    low_data = low_list{i};
    
    % Real data
    xu = up_data(:,1) * c0;
    yu = up_data(:,2) * c0;
    xl = low_data(:,1) * c0;
    yl = low_data(:,2) * c0;

    % Apply smoothing
    xu_smooth = movmean(up_data(:,1) * c0, smooth_win);
    yu_smooth = movmean(up_data(:,2) * c0, smooth_win);
    xl_smooth = movmean(low_data(:,1) * c0, smooth_win);
    yl_smooth = movmean(low_data(:,2) * c0, smooth_win);

    % Compute camber line, rotated so that the chord is horizontal
    % Chord length used in Prandtl's theory
    [~, ~, ~, ~, c_i] = compute_camber_MP_rotated(xu, yu, xl, yl, Npts);

    % In thin airfoi we use smooth data due to manual inprecision
    [M_i, P_i, xc_i, yc_i] = compute_camber_MP_rotated(xu_smooth, yu_smooth, xl_smooth, yl_smooth, Npts);

    % Store results
    xc_cells{i} = xc_i;
    yc_cells{i} = yc_i;
    M(i) = M_i;
    P(i) = P_i;
    c_values(i) = c_i;
end



%% Thin airfoil theory

alpha_plot_deg = linspace(-10, 10, 300);   % [deg]
alpha_plot_rad = deg2rad(alpha_plot_deg);  % [rad]

nCases = numel(M);

% Matrix initialization
A0_TAT = zeros(1, nCases);
A1_TAT = zeros(1, nCases);
A2_TAT = zeros(1, nCases);
alphaL0_2D = zeros(1, nCases);
Cm_c4_2D = zeros(1, nCases);
a_infty = zeros(1, nCases);
alpha_L0_thin = zeros(1, nCases);
Cl_TAT = zeros(nCases, numel(alpha_plot_rad));

for i = 1:nCases

    % Actual camber line
    xc_i = xc_cells{i};       % x/c
    yc_i = yc_cells{i};       % camber line

    % dyc/dx
    dyc_dx_vals = gradient(yc_i) ./ gradient(xc_i);
    dyc_dx_fun = @(x) interp1(xc_i*c0, dyc_dx_vals, x, 'linear', 'extrap');

    % Compute A0, A1, A2
    [A0, An] = thin_airfoil_coeffs(dyc_dx_fun, c0, 2, 2000);
    A0_TAT(i) = A0;
    A1_TAT(i) = An(1);
    A2_TAT(i) = An(2);

    % Cl(alpha), alpha_L0, Cm_c/4
    [Cl_i, alphaL0_2D(i), Cm_c4_2D(i)] = thin_airfoil_CL(alpha_plot_rad, A0, A1_TAT(i), A2_TAT(i));
    Cl_TAT(i,:) = Cl_i;

    % 2D slope
    p_fit = polyfit(alpha_plot_rad, Cl_i, 1);
    a_infty(i) = p_fit(1);
    alpha_L0_thin(i) = -p_fit(2)/p_fit(1);
end


%% Plot Thin Airfoil Theory vs Experimental Data

fprintf('\n=== Thin Airfoil Theory ===\n');
fprintf('Slope = %.3f\n [1/rad]', 2*pi)
figure; hold on; grid on; box on;

alpha_L0_thin_deg = rad2deg(alphaL0_2D);  

% Experimental data (discrete points)
h_exp = plot(alpha_fit_deg, CL_fit, 'ko','LineWidth', 1.2, 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Experimental (linear region)');

% Thin Airfoil Theory curves
colors = lines(nCases);
h_hat = gobjects(nCases,1);

for i = 1:nCases

    h_hat(i) = plot(alpha_plot_deg, Cl_TAT(i,:), 'LineWidth', 1.8, 'Color', colors(i,:));
    
    if i == 1
        set(h_hat(i), 'DisplayName', 'Flat plate');
        fprintf('alpha_L0 for standard case (without detachment) = %.3f [deg]\n', ...
                alpha_L0_thin_deg(i));
    else
        set(h_hat(i), 'DisplayName', sprintf('Profile at %.0f degree', 5*i));
        fprintf('alpha_L0 for image at %.0f degree = %.3f [deg]\n', 5*i, alpha_L0_thin_deg(i));
    end
end

xlabel('\alpha [deg]');
ylabel('C_l');
title('Thin Airfoil Theory vs Experimental Fit');

xline(0, 'k-', 'LineWidth', 1);
yline(0, 'k-', 'LineWidth', 1);

axis tight;

legend([h_exp; h_hat(:)], 'Location', 'NorthWest');

%saveas(gcf, 'ThinAirfoilTh.jpg')


%% Prandtl: ideal elliptical case: delta = 0
a_infty = a_infty(1);       % 2pi from thin airfoil theory
%a_infty = 4.6;              % Value extrapolated from Figliozzi et al. 2018 – "Lift curve slope variation at low Reynolds numbers" 
delta_ell = 0;
a_ell = a_infty / ( 1 + (a_infty/(pi*AR)) * (1 + delta_ell) );  % [1/rad]

%% General case: obtain delta_eff by setting a = a_exp
delta_eff = ( (a_infty/a_exp - 1) * (pi*AR/a_infty) ) - 1;


%% Geometric delta from the lifting line
a_2d = a_infty;         % 2D slope used in LL [1/rad]
alpha_geo_deg = 8;      % Angle of attack (arbitrary choice)
alpha_geo_rad = deg2rad(alpha_geo_deg);

delta_geo = zeros(6,1);
a_geo = zeros(6,1);

for i = 1:6
    % Chord and twist functions
    c_fun = @(y) c_values(i);
    twist_fun = @(y) 0;

    % Compute delta from geometry
    [delta_geo(i), A] = compute_delta_from_geometry(b, c_fun, twist_fun, alpha_geo_rad, a_2d, alpha_L0_thin(i));

    a_geo(i) = a_infty / ( 1 + (a_infty/(pi*AR)) * (1 + delta_geo(i)) );
end 

%% Results comparison
% Plot the results obtained to obtain a numerical comparison
fprintf('\n=== Slope comparison ===\n');
fprintf('Slope measured in the paper = %.3f  [1/rad]\n', a_exp);
fprintf('Delta from experimental data = %.3f\n', delta_eff);
fprintf('Slope thin airfoil theory = %.3f  [1/rad]\n', 2*pi);
fprintf('Slope ellittic (delta=0) = %.3f  [1/rad]\n', a_ell);

angles = [0, 10, 15, 20, 25, 30];
for i = 1:6
    fprintf('\n--- Geometry extracted from %d° image ---\n', angles(i));
    fprintf('Delta geometric (%d°) = %.4f\n', angles(i), delta_geo(i));
    fprintf('Slope geometric (%d°) = %.4f [1/rad]\n', angles(i), a_geo(i));
end



%% Delta taken from Glauert 1959

 delta_gla = 0.12;
 a_gla = a_infty / ( 1 + (a_infty/(pi*AR)) * (1 + delta_gla) );
 fprintf('\nDelta Glauert = %.3f\n', delta_gla);
 fprintf('Slope Glauert = %.3f [1/rad]\n', a_gla);


%% Plot of results

CL_ell = a_ell .* (alpha_plot_rad - alphaL0);
CL_delta = a_exp .* (alpha_plot_rad - alphaL0);
CL_gla = a_gla .* (alpha_plot_rad - alphaL0);
CL_thin_airfoil = 2*pi .* (alpha_plot_rad - alphaL0);
CL_geo = zeros(6, numel(alpha_plot_rad));
for i = 1:6
    CL_geo(i,:) = a_geo(i) .* (alpha_plot_rad - alphaL0);
end

figure; hold on; grid on; box on;

% Theoretical curves
h1 = plot(alpha_plot_deg, CL_ell, '-', 'LineWidth', 1.8);
h2 = plot(alpha_plot_deg, CL_geo(1,:), '-.', 'LineWidth', 1.8);
h3 = plot(alpha_plot_deg, CL_thin_airfoil, '-.', 'LineWidth', 1.8);

% Experimental data
h4 = plot(alpha_fit_deg, CL_fit, 'o', 'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', 'auto');

% Cartesian axes
xline(0, 'k-', 'LineWidth', 1);
yline(0, 'k-', 'LineWidth', 1);

xlabel('\alpha [deg]');
ylabel('C_L');

if a_infty == 4.6 
    title('C_L(\alpha): Elliptical vs \delta_{geo} vs Thin Airfoil for a_{2D} = 4.6 [1/rad]');
else
    title('C_L(\alpha): Elliptical vs \delta_{geo} vs Thin Airfoil for a_{2D} = 2\pi [1/rad]');
end

legend_text_geo = sprintf('Geometric (\\delta_{geo}=%.3f)', delta_geo(1));

legend([h1 h2 h3 h4], {'Elliptical (\delta=0)', legend_text_geo, 'Thin airfoil', 'Experimental'}, ...
       'Location', 'NorthWest');

axis tight;

%saveas(gcf, 'PrandtlTh46.jpg')


%% ---------------- AUXILIARY FUNCTION ----------------
function [delta, A_n] = compute_delta_from_geometry(b, c_fun, twist_fun, alpha_root, a_2d, alpha0)
% Function that compute delta from the geomety based on Prandtl's theory

% INPUT
% b: total wingspan [m]
% c_fun(y): local chord [m]
% twist_fun(y): local geometric twist [rad]
% alpha_root: root angle [rad]
% a_2d: slope 2D [1/rad]

% OUTPUT
% delta: deviation from the elliptical profile [-]
% A_n: matrix containing the coefficient of A_n
    
    N = 501;                               % number of stations on the semi-opening
    theta = linspace(0.001, pi-0.001, N);  % Prandtl angle
    y = -(b/2) * cos(theta);               % spanwise coordinates on half wing
    
    % Local distributions
    c = arrayfun(c_fun, abs(y));
    twist = arrayfun(twist_fun, abs(y));
    alpha_local = alpha_root + twist - alpha0;
    
    N_modes = 501;

    % Prandtl system matrix
    A_matrix = zeros(N, N_modes);
    % Right hand side of the equation
    RHS = zeros(N, 1);
    
    for i = 1:N
        for n = 1:N_modes

            A_matrix(i,n) = sin(n*theta(i)) * (1+ c(i)*a_2d*n/(4*b*sin(theta(i))));
    
        end
        RHS(i) = alpha_local(i)*c(i)*a_2d/(4*b);
    end
    
    % I solve for the coefficients A_n
    A_n = A_matrix \ RHS;
    
    % Delta calculation
    delta = 0;
    for n = 2:N_modes
        delta = delta + n * (A_n(n)/A_n(1))^2;
    end
end

function [A0, An] = thin_airfoil_coeffs(dyc_dx_fun, c, N_modes, N_theta)
% Compute A0 and An from the thin profile theory using numerical integration.

% INPUT
% dyc_dx_fun: handle
% c : Chord [m]
% N_modes: Number of coefficients An
% N_theta: Number of integration points on theta in [0,pi]

% OUTPUT
% A0: Scalar
% An: Vector [N_modes x 1]

    % mapping x in function of theta
    theta = linspace(0, pi, N_theta);
    x = 0.5 * c * (1 - cos(theta));
    
    dyc_dx = dyc_dx_fun(x);

    % A0
    A0 = (1/pi) * trapz(theta, dyc_dx);

    % An
    An = zeros(N_modes,1);
    for n = 1:N_modes
        An(n) = (2/pi) * trapz(theta, dyc_dx .* cos(n*theta));
    end
end

function [Cl, alpha_L0, Cm_c4] = thin_airfoil_CL(alpha_rad, A0, A1, A2)
% Returns Cl, alpha_L0, and Cm_c/4 from the thin airfoil theory.

% INPUT
% alpha_rad: angle of attack [rad]
% A0, A1, A2: camber line coefficients

% OUTPUT
% Cl: 2D lift coefficient
% alpha_L0: zero lift angle [rad]
% Cm_c4: Moment coefficient at 1/4 of the string length

    alpha_L0 = (2*A0 - A1)/2;
    Cl = 2*pi*alpha_rad + pi*(A1 - 2*A0);
    Cm_c4 = -pi/4 * (A1 - A2);
end

function [M, P, xc_grid, yc, c] = compute_camber_MP(xu, yu, xl, yl, Npts)
% Calculate M (maximum camber) and P (position of maximum camber) from the 
% coordinates of the upper and lower surfaces of a wing profile.

% INPUT
% xu, yu: upper surface coordinates
% xl, yl: lower surface coordinates
% Npts: number of points along the string for discretization

% OUTPUT
% M: maximum camber [-]
% P: x/c position of maximum camber [-]
% xc_grid: dimensionless x/c grid
% yc: dimensionless camber line
% c: chord length

    % Put everything in a vector column
    xu = xu(:); yu = yu(:);
    xl = xl(:); yl = yl(:);

    % Order the points based on x
    [xu, idxu] = sort(xu);
    yu = yu(idxu);

    [xl, idxl] = sort(xl);
    yl = yl(idxl);

    % LE, TE and the chord
    x_all = [xu; xl];
    x_le = min(x_all);
    x_te = max(x_all);
    c = x_te - x_le;

    % Uniform grid along the chord
    xc_grid = linspace(0, 1, Npts)';      % x/c adimensional
    x_phys = x_le + c * xc_grid;          % physical coordinates x

    % Interpolate the top and bottom surfaces on the same grid
    yu_i = interp1(xu, yu, x_phys, 'linear', 'extrap');
    yl_i = interp1(xl, yl, x_phys, 'linear', 'extrap');

    % Compute the camber line
    yc_phys = 0.5 * (yu_i + yl_i);
    yc = yc_phys / c;

    % Find the maximum camber and its position
    [M, idx_max] = max(yc);     % M = y_c,max / c
    P = xc_grid(idx_max);       % P = (x/c)
end

function [M, P, xc_grid, yc, c] = compute_camber_MP_rotated(xu, yu, xl, yl, Npts)
% Place the profile at an angle so that the chord is horizontal.

% INPUT
% xu, yu: upper surface coordinates (rotated profile)
% xl, yl: lower surface coordinates (rotated profile)
% Npts: number of points along the chord for discretization

% OUTPUT
% M: maximum camber [-]
% P: x/c position of maximum camber [-]
% xc_grid: dimensionless x/c grid
% yc: dimensionless camber line
   
    % Column vectors
    xu = xu(:); yu = yu(:);
    xl = xl(:); yl = yl(:);

    % Determine LE and TE
    x_all = [xu; xl];
    y_all = [yu; yl];

    [~, idx_le] = min(x_all);
    x_le = x_all(idx_le);
    y_le = y_all(idx_le);

    [~, idx_te] = max(x_all);
    x_te = x_all(idx_te);
    y_te = y_all(idx_te);

    % Chord length
    c = x_te - x_le;

    % Rotation angle
    alpha = atan2(y_te - y_le, x_te - x_le);

    % Rotate all points
    xu_rot = (xu - x_le)*cos(-alpha) - (yu - y_le)*sin(-alpha);
    yu_rot = (xu - x_le)*sin(-alpha) + (yu - y_le)*cos(-alpha);

    xl_rot = (xl - x_le)*cos(-alpha) - (yl - y_le)*sin(-alpha);
    yl_rot = (xl - x_le)*sin(-alpha) + (yl - y_le)*cos(-alpha);

    % Compute camber etc.
    [M, P, xc_grid, yc] = compute_camber_MP(xu_rot, yu_rot, xl_rot, yl_rot, Npts);
end
