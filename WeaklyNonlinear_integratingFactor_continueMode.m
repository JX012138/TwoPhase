clearvars; close all; clc;
tic;

%% 0. RESTART CONFIGURATION
restart_mode = true;       % Set to TRUE to continue from a previous file
restart_time = 15.0;       % The time of the file you want to load (e.g., 30.0)
T_final      = 16.0;       % The new final time you want to reach

%% 1. Parameters & Precomputation
% parameters 1
m = 0.2; 
Re = 6; 
W_rod = 0.0; 
wavenumber_target = 0.2;
p_normalmode = 1;
L = p_normalmode*pi/wavenumber_target;
Ca = 10;

% parameters 2
a_bar = 0.1;  
delta = 0.1; 
save_interval = 0.01; 
nu = (pi/L)^2;

% --- Base-flow Constants ---
A1 = (1 - a_bar^2 - 4*W_rod) / (4*log(1/a_bar));
C = (m-1)/(2*m) + (1-m)/m * A1; 
G2 = 1 / (3*m*Ca); 

% --- Spectral Resolution (Half-Spectrum) ---
Nz = 256; 
Nz_half = Nz/2 + 1; 
p_modes_half = 0:(Nz/2); 
dealias_mask_half = ones(size(p_modes_half));
p_max = Nz/3; 
dealias_mask_half(p_modes_half > p_max) = 0;

% --- Precompute Linear Operator (Same as before) ---
N_SP = 66; 
[D1, x_cheb] = cheb(N_SP);
alpha = (1 - a_bar)/2;
beta = (1 + a_bar)/2;
r = alpha*x_cheb + beta;
B1 = W_rod - A1*log(a_bar) + 0.25*a_bar^2;
w0 = -0.25*r.^2 + A1*log(r) + B1;
w0_r = -1/2 *r + A1./r;
w0_rr = -1/2 - A1./(r.^2);
G3_modes_half = zeros(size(p_modes_half));
for idx = 1:length(p_modes_half)
    p = p_modes_half(idx);
    p_new = (pi/L)*p; 
    zeta_rr_at_1 = compute_zeta_rr(p_new, Re, N_SP, D1, r, alpha, w0, w0_r, w0_rr);
    G3_modes_half(idx) = (1i * p_new) * C * (1 - zeta_rr_at_1) / (2 * m); 
end
lin_coeff_half = (p_modes_half.^2 - nu * (p_modes_half.^4)) + G3_modes_half/(G2*nu);  

%% 2. Time-step (Integrating Factor)
dt_main = 5e-4; % Precision step
E = exp(lin_coeff_half * dt_main);
E_half = exp(lin_coeff_half * dt_main/2);

%% 3. INITIAL CONDITION (From Current Folder)
x = linspace(-pi,pi,Nz+1); 
x = x(1:end-1);

if restart_mode
    % --- RESTART LOGIC ---
    % Look for file in the CURRENT directory
    fname = sprintf('newrescaled_h_at_t%.3f.dat', restart_time);
    
    if ~exist(fname, 'file')
        error('Restart file NOT found in current folder: %s', fname);
    end
    
    fprintf('RESTARTING: Loading %s ...\n', fname);
    data = readmatrix(fname);
    h = data(:, 2)'; % Extract h column and transpose to row
    
    % Validation
    if length(h) ~= Nz
        error('Grid mismatch! File has %d points, Code expects %d.', length(h), Nz);
    end
    
    current_start_time = restart_time;
else
    % --- FRESH START ---
    fprintf('FRESH START: Initializing t=0...\n');
    h = delta * cos(p_normalmode*x); 
    current_start_time = 0.0;
    
    % Save Initial Condition if fresh start
    filename = 'newrescaled_h_at_t0.000.dat';
    writematrix([x', h'], filename);
    fprintf('Saved %s\n', filename);
end

% FFT Init
h_hat_full = fft(h);
h_hat_half = h_hat_full(1:Nz_half);

%% 4. Production Run
% Calculate steps needed
if T_final <= current_start_time
    fprintf('Simulation already reached target time %.1f. Exiting.\n', T_final);
    return;
end

Nt = round((T_final - current_start_time) / dt_main);
fprintf('Running %d steps (t=%.3f -> %.3f)\n', Nt, current_start_time, T_final);

next_save_t = current_start_time + save_interval;
Coeff_NL = 1.0; 

for n = 1:Nt
    % IFRK4 Step
    k1 = compute_nonlinear_half(h_hat_half, p_modes_half, dealias_mask_half, Coeff_NL);
    temp = (h_hat_half + (dt_main/2) * k1) .* E_half;
    k2 = compute_nonlinear_half(temp, p_modes_half, dealias_mask_half, Coeff_NL);
    temp = h_hat_half .* E_half + (dt_main/2) * k2;
    k3 = compute_nonlinear_half(temp, p_modes_half, dealias_mask_half, Coeff_NL);
    temp = h_hat_half .* E + dt_main * k3 .* E_half;
    k4 = compute_nonlinear_half(temp, p_modes_half, dealias_mask_half, Coeff_NL);
    
    h_hat_half = h_hat_half .* E + (dt_main/6) * ...
                 (k1 .* E + 2*k2 .* E_half + 2*k3 .* E_half + k4);
    h_hat_half = h_hat_half .* dealias_mask_half;
    
    % --- Current Time Calculation ---
    current_t = current_start_time + n * dt_main;
    
    % --- Save Logic (Directly to Current Folder) ---
    if current_t >= next_save_t - dt_main/2
        h_full = reconstruct_real_h(h_hat_half);
        if any(isnan(h_full)) || any(isinf(h_full)), warning('NaN/Inf at step %d', n); break; end

        % Save simply as 'rescaled_h_at_tX.XXX.dat'
        filename = sprintf('newrescaled_h_at_t%.3f.dat', next_save_t);
        writematrix([x', h_full'], filename);
        fprintf('Saved %s\n', filename);
        
        next_save_t = next_save_t + save_interval;
        if next_save_t > T_final, break; end
    end
end
disp('Simulation complete.');
toc;

%% --- Helper Functions (Standard) ---
function nk_half = compute_nonlinear_half(h_hat_half, p_half, mask_half, Coeff)
    h_hat_full = [h_hat_half, conj(h_hat_half(end-1:-1:2))];
    h_phys = real(ifft(h_hat_full));
    p_full = [p_half, -p_half(end-1:-1:2)];
    h_x_phys = real(ifft(1i * p_full .* h_hat_full));
    nl_phys = h_phys .* h_x_phys;
    nl_hat_full = fft(nl_phys);
    nk_half = -Coeff * nl_hat_full(1:length(p_half)) .* mask_half;
end

function h = reconstruct_real_h(h_hat_half)
    h_hat_full = [h_hat_half, conj(h_hat_half(end-1:-1:2))];
    h = real(ifft(h_hat_full));
end

function [D, x] = cheb(N)
    if N == 0, D = 0; x = 1; return; end
    x = cos(pi * (0:N)' / N);
    c = [2; ones(N-1,1); 2] .* (-1).^(0:N)';
    X = repmat(x, 1, N+1);
    dX = X - X';
    D = (c * (1./c)') ./ (dX + eye(N+1));
    D = D - diag(sum(D, 2));
end

function zrr = compute_zeta_rr(k_phys, Re, N_SP, D1, r, alpha, w0, w0_r, w0_rr)
    M0 = diag(k_phys^4 + 1i*k_phys*Re*w0_rr - (1i*k_phys*Re*w0_r)./r + 1i*k_phys^3*Re*w0);
    M1 = diag((-3./r.^3 + 2*k_phys^2./r + 1i*k_phys*Re*w0./r) / alpha);
    M2 = diag((3./r.^2 - 2*k_phys^2 - 1i*k_phys*Re*w0) / alpha^2);
    M3 = diag(-2 ./ (r * alpha^3));
    M4 = diag(ones(N_SP+1, 1) / alpha^4);
    A = M4 * D1^4 + M3 * D1^3 + M2 * D1^2 + M1 * D1 + M0;
    A(1,:) = 0; A(1,1) = 1; A(2,:) = D1(1,:);
    A(end,:) = 0; A(end,end) = 1; A(end-1,:) = D1(end,:);
    b = zeros(N_SP+1, 1); b(2) = alpha;
    zeta = A \ b;
    zeta_rr_xx = D1^2 * zeta;
    zrr = zeta_rr_xx(1) / (alpha^2);
end