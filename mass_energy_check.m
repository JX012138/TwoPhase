clearvars; close all; clc;

%% 1. Find and sort files
files = dir('rescaled_h_at_t*.dat');
if isempty(files)
    error('No files found matching rescaled_h_at_t*.dat');
end

% Extract times from filenames
time_points = zeros(1, length(files));
for i = 1:length(files)
    filename = files(i).name;
    
    % Robust time extraction
    t_start = strfind(filename, 't');
    if length(t_start) >= 2
        % Format: rescaled_h_at_tX.XXX.dat
        time_str = filename(t_start(2)+1:end-4);  % Remove '.dat'
    else
        % Fallback: try regex
        time_match = regexp(filename, 't(\d+\.?\d*)\.dat', 'tokens');
        if ~isempty(time_match)
            time_str = time_match{1}{1};
        else
            error('Cannot parse time from: %s', filename);
        end
    end
    time_points(i) = str2double(time_str);
end

% Sort by time
[time_points, sort_idx] = sort(time_points);
files = files(sort_idx);

masses = zeros(1, length(files));
h_mean = zeros(1, length(files));  
energies = zeros(1, length(files));

for i = 1:length(files)
    filename = files(i).name;
    data = readmatrix(filename);
    
    if size(data, 2) == 2
        x = data(:, 1);
        h = data(:, 2);
    else
        error('File %s has unexpected format', filename);
    end
    
    dx = x(2)-x(1);

    % Check for NaN/Inf
    if any(isnan(h)) || any(isinf(h))
        warning('File %s contains NaN/Inf values!', filename);
    end
    
    x = [x; x(end)+dx];
    h = [h; h(1)];

    % Compute mass using trapezoidal rule
    masses(i) = trapz(x, h);
    energies(i) = trapz(x, h.^2);

end

%% 3. Plotting
figure('Position', [100, 100, 1000, 800]);

% Subplot 1: Mass (Should be CONSTANT)
subplot(2,1,1);
% Plot relative error to see machine precision
plot(time_points, masses, 'b-', 'LineWidth', 2);
title('Mass Conservation Check');
xlabel('Time'); ylabel('Mass');
grid on;
subtitle('Ideally should be close to 0 ');

% Subplot 2: Energy (Should GROW then SATURATE)
subplot(2,1,2);
plot(time_points, energies, 'b-', 'LineWidth', 2);
title('System Energy (L_2 Norm)');
xlabel('Time'); ylabel('Energy E(t) = \int h^2 dx');
grid on;
subtitle('Expectation: Exponential growth followed by saturation');