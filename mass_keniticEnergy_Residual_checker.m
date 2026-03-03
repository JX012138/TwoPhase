clearvars; close all; clc;

%% 1. Parameters (Must match WeaklyNonlinear_integratingFactor_continueMode.m)
m = 0.2; 
Re = 6; 
W_rod = 0.0; 
wavenumber_target = 0.2; 
p_normalmode = 1;
L = p_normalmode*pi/wavenumber_target; 
Ca = 10; 
a_bar = 0.1;
nu = (pi/L)^2;
A1 = (1 - a_bar^2 - 4*W_rod) / (4*log(1/a_bar));
C = (m-1)/(2*m) + (1-m)/m * A1; 
G2 = 1 / (3*m*Ca); 

%% 2. File Handling
files = dir('newrescaled_h_at_t*.dat');
if isempty(files), error('No data files found.'); end

time_points = zeros(1, length(files));
for i = 1:length(files)
    t_start = strfind(files(i).name, 't');
    time_points(i) = str2double(extractBetween(files(i).name, t_start(2)+1, '.dat'));
end
[time_points, sort_idx] = sort(time_points);
files = files(sort_idx);

%% 3. Precompute Linear Operator G3 (Spectral)
Nz = 256; % Must match simulation resolution
p = [0:Nz/2-1, 0, -Nz/2+1:-1]'; % Wavenumbers

%% 4. Diagnostic Loop
mass_drift = zeros(size(files));
L2_norm = zeros(size(files));
convective_residual = zeros(size(files));
energy_rhs = zeros(size(files));

for i = 1:length(files)
    data = readmatrix(files(i).name);
    x = data(:, 1); h = data(:, 2);
    dx = x(2) - x(1);
    
    % --- Mass Conservation ---
    mass_drift(i) = trapz([x; x(end)+dx], [h; h(1)]);
    
    % --- L2 Norm (Energy) ---
    L2_norm(i) = 0.5 * trapz([x; x(end)+dx], [h; h(1)].^2);
    
    % --- Convective Residual Check ---
    % Analytical: \int h^2 * h_x dx = 0
    h_hat = fft(h);
    hx = real(ifft(1i * p .* h_hat));
    convective_residual(i) = trapz([x; x(end)+dx], ( [h; h(1)].^2 .* [hx; hx(1)]));
    
    % --- RHS Energy Budget ---
    % Term 1: hx^2 (Surface Tension / Growth)
    % Term 2: nu*hxx^2 (Dissipation)
    hxx = real(ifft(-(p.^2) .* h_hat));
    term_linear = trapz([x; x(end)+dx],[hx; hx(1)].^2 - nu*[hxx;hxx(1)].^2);
    
    % Term 3: Non-local G3 term (Linearized Hydrodynamics)
    % Calculation: \int h * \mathcal{F}^{-1}(G3 * h_hat) dx
    % (Implementation requires G3_vector from solver)
    
    energy_rhs(i) = term_linear; % Update with G3 and Scaling factors
end

%% 5. Visualization
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1); plot(time_points, mass_drift, 'k-'); title('Mass Drift');
subplot(1,3,2); plot(time_points, convective_residual, 'r-'); title('Convective Residual');
subplot(1,3,3); plot(time_points, L2_norm, 'b-'); title('Energy E(T)');