% MATLAB script to plot max and min vertical positions from .dat files
clc;clear;
% Get list of files
files = dir('rescaled_h_at_t*.dat');
num_files = length(files);

%parameters
N = 66;
a_bar = 0.1;
m = 2; 
Re = 1000; 
W_rod = -0.1; 
Ca = 10;
L = 2.094;
nu = (pi/L)^2;
p_normalmode = 1;
p = p_normalmode;
p_sqrtnu = p * sqrt(nu);
delta = 0.1; 
wavenumber = p_normalmode*pi/L;
fprintf("wavenumber = %.2f \n", wavenumber);

% Base-flow constants
A1 = (1 - a_bar^2 - 4*W_rod) / (4*log(1/a_bar));
B1 = W_rod - A1 * log(a_bar) + (1/4) * a_bar^2;
C = (m-1)/(2*m) + (1-m)/m * A1; % bracket term in G3

% Scaling coefficients
G1 = (A1/m - 1/(2*m));
G2 = 1 / (3*m*Ca); % prefactor in nonlocal term

% Preallocate arrays
times_compare = zeros(num_files, 1);
max_vals = zeros(num_files, 1);
min_vals = zeros(num_files, 1);

% Loop over files
for i = 1:num_files
    filename = files(i).name;
    
    % Extract time from filename (string between the second 't' and '.dat')
    % Find positions of all 't' in the filename
    t_indices = strfind(filename, 't');
    % Use the second 't' as the start point
    time_str = extractBetween(filename, t_indices(2)+1, '.dat');
    T = str2double(time_str);
    times_compare(i) = T;

    % Load data from file
    data = load(filename);
    
    if isempty(data) || size(data, 2) ~= 2
        fprintf('Warning: %s is empty or has incorrect format. Skipping.\n', filename);
        continue;
    end

    % Assume second column is the vertical coordinate
    y = data(:, 2);
    
    % Compute max and min
    max_vals(i) = max(y);
    min_vals(i) = min(y);
end

% Sort by time
[times_compare, idx] = sort(times_compare);
max_vals = max_vals(idx);
min_vals = min_vals(idx);

%% linear growth rate

% Define constants from the transformation r = alpha x + beta
alpha = (1 - a_bar) / 2;
beta = (1 + a_bar) / 2;

% Compute differentiation matrices (this will give nodes x from 1 to -1
[D1, x] = cheb(N);
D2 = D1 * D1;
D3 = D1 * D2;
D4 = D1 * D3;

% Initialize coefficient matrices
M0 = zeros(N+1);
M1 = zeros(N+1);
M2 = zeros(N+1);
M3 = zeros(N+1);
M4 = zeros(N+1);

% Compute w0(x), w0_x(x), w0_xx(x)
w0 = - (1/4) * (alpha * x + beta).^2 + A1 * log(alpha * x + beta) + B1;
w0_x = (-1/2) * (alpha * x + beta)  + (A1 ) ./ (alpha * x + beta);
w0_xx = - 1 / 2 - (A1 ) ./ (alpha * x + beta).^2;

% Evaluate coefficients at each node
for j = 1:N+1
    
    % Collect the corresponding base velocity and its derivatives 
    w0_j = w0(j);
    w0_x_j = w0_x(j);
    w0_xx_j = w0_xx(j);

    % Compute coefficient matrices Mn(j,j) for different order
    M0(j,j) = p_sqrtnu^4 + 1i * p_sqrtnu * Re * w0_xx_j - (1i *p_sqrtnu * Re * w0_x_j) / (alpha * x(j) + beta) + 1i * p_sqrtnu^3 * Re * w0_j;
    M1(j,j) =  (-3 / (alpha * x(j) + beta)^3 + (2 * p_sqrtnu^2) / (alpha * x(j) + beta) + (1i * p_sqrtnu * Re * w0_j) / (alpha * x(j) + beta)) * (1 / alpha); 
    M2(j,j) = (3 / (alpha * x(j) + beta)^2 - 2 * p_sqrtnu^2 - 1i * p_sqrtnu * Re * w0_j ) * (1 / alpha^2);
    M3(j,j) = - ( 2 / (alpha * x(j) + beta) ) * (1 / alpha^3);
    M4(j,j) = (1 / alpha^4);
end

% Initialize A and b
A = M4 * D4 + M3 * D3 + M2 * D2 + M1 * D1 + M0;
b = zeros(N+1, 1);

% Apply boundary conditions by modifying A and b
% At x = x(1) = 1
A(1, :) = 0;
A(1, 1) = 1;  % Enforce zeta(1) = 0
b(1) = 0;

A(2, :) = D1(1, :);  % Enforce zeta_x(1) = alpha
b(2) = alpha;

% At x = x(N+1) = -1
A(N+1, :) = 0;
A(N+1, N+1) = 1;  % Enforce zeta(N+1) = 0
b(N+1) = 0;

A(N, :) = D1(N+1, :);  % Enforce zeta_x(-1) = 0
b(N) = 0;

% Solve the linear system
zeta = A \ b;

% Compute the second derivative
zeta_dd = D2 * zeta;

% Evaluate at x = 1 (first node)
zeta_dd_at_1 = zeta_dd(1);

% transform the coordinate back to r=1
zetar_dd_at_1 = (1/alpha)^2 * zeta_dd_at_1;

% Extract the imaginary part
Im_zeta_dd_at_1 = imag(zetar_dd_at_1);

% input the function of lamda and taking the real part of it
Re_lamda = (p^2 - nu*p^4) + (1/(G2*nu))*((p*sqrt(nu))/(2*m)) *C*Im_zeta_dd_at_1;
growth_rate = delta*exp(Re_lamda.*times_compare);


% Convert to log scale (you're already doing this for plotting)
log_amplitudes = log(max_vals);

% Find saturation amplitude and time
A_sat = exp(mean(log_amplitudes(end-10:end))); % Convert back from log

% Find saturation time (when growth effectively stops)
% Method: When the slope becomes very small
slopes = diff(log_amplitudes) ./ diff(times_compare(1:end));
saturation_idx = find(abs(slopes) < 0.1 * max(abs(slopes)), 1);  %%%%% the limit of this saturation point can be changed
t_sat = times_compare(saturation_idx);

fprintf('Saturation Analysis:\n');
fprintf('Saturation Amplitude: %.6f\n', A_sat);
fprintf('Saturation Time: %.3f\n', t_sat);


% Plot
figure;
plot(times_compare, log(max_vals), 'b-', 'LineWidth', 1.5);
hold on;
plot(times_compare, log(growth_rate), 'r-', 'LineWidth', 1.5);
plot(t_sat, log_amplitudes(saturation_idx), 'o', 'LineWidth',1.5);
hold off;
xlabel('Time');
ylabel('Position');
title('Max Position and Growth Rate vs Time');
legend('Max Position', 'Growth Rate', 'Saturation Point');
grid on;


