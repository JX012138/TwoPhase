clearvars; close all; clc;

verticalOffset = 8;                   % Adjust for spacing between interfaces

files = dir('rescaled_h_at_t*.dat');

% Sort files by time
time_points = zeros(1, length(files));
for i = 1:length(files)
    filename = files(i).name;
    t_indices = strfind(filename, 't');
    time_str = extractBetween(filename, t_indices(2)+1, '.dat');
    time_points(i) = str2double(time_str);
end

[sorted_times, sort_idx] = sort(time_points);
sorted_files = files(sort_idx);

user_input = input('Enter desired times (e.g., [0.4, 0.8, 1.0]): ');
desired_times = user_input;

% Find closest matches
selected_indices = [];
selected_times = [];
tolerance = 0.05;

for i = 1:length(desired_times)
    [min_diff, idx] = min(abs(sorted_times - desired_times(i)));
    if min_diff <= tolerance
        selected_indices = [selected_indices, idx];
        selected_times = [selected_times, sorted_times(idx)];
    end
end
    
% Remove duplicates and sort
[selected_indices, unique_idx] = unique(selected_indices, 'stable');
selected_times = selected_times(unique_idx);
[selected_times, sort_order] = sort(selected_times);
selected_indices = selected_indices(sort_order);

% Get selected files
selectedFiles = sorted_files(selected_indices);

fprintf('Found %d .dat files\n', length(sorted_files));
fprintf('Selected %d files for plotting:\n', length(selectedFiles));
for i = 1:length(selectedFiles)
    fprintf('  inte_%d: %s\n', i, selectedFiles(i).name);
end

% Create figure
figure('Position', [100, 100, 1200, 800]);

% Color scheme
colors = parula(length(selectedFiles));  % or use: colors = jet(length(selectedFiles));

% Pre-read first file to get x coordinates
firstData = load(fullfile(selectedFiles(1).name));
if size(firstData, 2) >= 2
    x = firstData(:, 1);  % First column is x
else
    x = (1:size(firstData, 1))';  % Use index if single column
end

hold on;

% Plot each selected interface
for i = 1:length(selectedFiles)
    filePath = fullfile(selectedFiles(i).name);
    data = load(filePath);
    
    % Extract height data
    if size(data, 2) >= 2
        h = data(:, 2);  % Second column is height
    else
        h = data(:, 1);  % Single column data
    end
    
    % Apply vertical offset for stacking
    h_offset = h + ((i-1) * verticalOffset);
    
    % Plot interface
    plot(x, h_offset, ...
         'LineWidth', 1.5, ...
         'Color', colors(i, :), ...
         'DisplayName', sprintf('inte_%d (t=%d)', i, time_points(sort_idx(selected_indices(i)))));
    
    % Add horizontal separator line
    plot([min(x), max(x)], [(i-1)*verticalOffset, (i-1)*verticalOffset], ...
         'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Customize plot
xlabel('Domain (x)', 'FontSize', 12);
ylabel('Interface Height (with offset)', 'FontSize', 12);
title(sprintf('Interface Evolution'), 'FontSize', 14);

% Create custom y-axis ticks for interface labels
yticks((0:length(selectedFiles)-1) * verticalOffset);
yticklabels(arrayfun(@(i) sprintf('inte_%d', i), 1:length(selectedFiles), 'UniformOutput', false));

% Add grid
grid on;
grid minor;
set(gca, 'GridAlpha', 0.3);

% Add legend
legend('Location', 'eastoutside', 'FontSize', 9);

% Adjust layout
axis tight;
box on;