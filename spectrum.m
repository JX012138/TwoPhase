clc;
clear;

files = dir('rescaled_h_at_t*.dat');
plot_n = 5;

% Extract and sort the first set of files by time
time_points1 = zeros(1, length(files));
for i = 1:length(files)
    filename = files(i).name;
    % Extract time from filename 
    t_indices = strfind(filename, 't');
    time_str = extractBetween(filename, t_indices(2)+1, '.dat');
    time_points1(i) = str2double(time_str);
end
[sorted_times1, sort_idx1] = sort(time_points1);
sorted_files1 = files(sort_idx1);

% Use the shorter length to avoid index issues
num_files = length(sorted_files1);

count = 1;
for n = num_files:-round(num_files/plot_n):1
    % Loop through each file in sorted order

    % Construct full filename from sorted lists
    filename = sorted_files1(n).name;
    
    % Read the data (two-column format: x y)
    data1 = readmatrix(filename);
    
    % Check if data is valid and has 2 columns
    if isempty(data1) || size(data1, 2) ~= 2
        fprintf('Warning: %s is empty or has incorrect format. Skipping.\n', filename);
        continue;
    end

    % Extract coordinates
    x1 = data1(:,1);  % Scale axial coordinate to [0, 2π] for MATLAB comparison
    y1 = data1(:,2);         % Radial coordinate (r)

    % Store time from sorted list
    times(count) = sorted_times1(n);

    % Get number of points
    Nx = length(x1);
    
    % Compute Fourier transform
    h_hat = fft(y1);
    
    % Generate wavenumber array (in correct order for FFT)
    k = [0:Nx/2-1, -Nx/2:-1]';  % wavenumbers
    k_positive = k(1:Nx/2);     % Only positive wavenumbers (0 to Nyquist)
    
    % Get amplitudes
    amplitude = abs(h_hat);
    amplitude_positive = amplitude(1:Nx/2);

    figure('Position', [100, 100, 800, 600]);
    semilogy(k_positive, amplitude_positive, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
    xlabel('Wavenumber $k$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$|\hat{h}(k)|$', 'Interpreter', 'latex', 'FontSize', 14);
    title(sprintf('Interface Perturbation Spectrum at $t = %.2f$', times(count)), ...
              'Interpreter', 'latex', 'FontSize', 16);
    grid on;
    xlim([0, max(k_positive)]);
    count = count +1;
end
