clearvars; close all; clc;

%% 1. Data Acquisition
files = dir('newrescaled_h_at_t*.dat');
time_points = zeros(length(files), 1);
for i = 1:length(files)
    t_idx = strfind(files(i).name, 't');
    time_points(i) = str2double(extractBetween(files(i).name, t_idx(2)+1, '.dat'));
end
[sorted_times, sort_idx] = sort(time_points);
sorted_files = files(sort_idx);

%% 2. Neighbor-Peak Tracking Logic (2D Selection)
n_snapshots = length(sorted_files); 
verticalOffset = 25; % Space between interfaces
selected_indices = zeros(n_snapshots, 1);
selected_indices(1) = find(sorted_times > 0, 1); % Start at saturation

current_idx = selected_indices(1);
for i = 2:n_snapshots
    % Identify current peak position
    data_curr = readmatrix(sorted_files(current_idx).name);
    [~, max_p] = max(data_curr(:,2));
    z_curr = data_curr(max_p, 1);
    
    found = false;
    % Search for next frame where peak has returned near z_curr
    for j = (current_idx + 10):length(sorted_files)
        data_next = readmatrix(sorted_files(j).name);
        [~, p_next] = max(data_next(:,2));
        
        % Selection criteria: Peak is spatially aligned with previous snapshot
        if abs(data_next(p_next, 1) - z_curr) < 0.05
            selected_indices(i) = j;
            current_idx = j;
            found = true;
            break;
        end
    end
    if ~found, selected_indices(i:end) = []; break; end
end

%% 3. 2D Stacked Plotting
figure('Position', [100, 100, 800, 900], 'Color', 'w');
hold on;
colors = winter(length(selected_indices)); % Gradient for time evolution

for i = 1:length(selected_indices)
    idx = selected_indices(i);
    if idx == 0, continue; end
    
    data = readmatrix(sorted_files(idx).name);
    z = data(:, 1);
    r = data(:, 2);
    
    % Apply vertical offset for stacking
    r_stacked = r + (i-1) * verticalOffset;
    
    % Plot interface
    plot(z, r_stacked, 'LineWidth', 1.8, 'Color', colors(i,:));
    
    % Add baseline for clarity (optional)
    plot([min(z), max(z)], [(i-1)*verticalOffset, (i-1)*verticalOffset], ...
         'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

%% 4. Aesthetics
xlabel('Spatial Domain $z$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Interface Height (Stacked)', 'Interpreter', 'latex', 'FontSize', 14);
title('2D Peak-Aligned Interface Evolution', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
axis tight;
set(gca, 'YTick', (0:n_snapshots-1)*verticalOffset);
set(gca, 'YTickLabel', arrayfun(@(t) sprintf('t=%.1f', t), ...
    sorted_times(selected_indices(selected_indices>0)), 'UniformOutput', false));