clearvars; close all; clc;

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

outputEveryN = 20;                    % Plot every 10th time step
verticalOffset = 8;                   % Adjust for spacing between interfaces

% Select every N-th file
selectedIdx = 1:outputEveryN:length(sorted_files);
selectedFiles = sorted_files(selectedIdx);

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
         'DisplayName', sprintf('inte_%d (t=%d)', i, time_points(sort_idx(selectedIdx(i)))));
    
    % Add horizontal separator line
    plot([min(x), max(x)], [(i-1)*verticalOffset, (i-1)*verticalOffset], ...
         'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Customize plot
xlabel('Domain (x)', 'FontSize', 12);
ylabel('Interface Height (with offset)', 'FontSize', 12);
title(sprintf('Interface Evolution (every %d time steps)', outputEveryN), 'FontSize', 14);

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
