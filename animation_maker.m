clearvars; close all; clc;

%% Simple interface animation
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

% Set up figure
figure('Position', [200, 200, 1000, 600]);

% Precompute global limits
all_data = cell(length(sorted_files), 1);
for i = 1:length(sorted_files)
    data = readmatrix(sorted_files(i).name);
    all_data{i} = data;
end

min_h = min(cellfun(@(d) min(d(:,2)), all_data));
max_h = max(cellfun(@(d) max(d(:,2)), all_data));
margin = 0.1 * (max_h - min_h);

% Animation loop
for i = 1:length(sorted_files)
    data = readmatrix(sorted_files(i).name);
    x = data(:,1);
    h = data(:,2);
    current_time = sorted_times(i);
    
    plot(x, h, 'b-', 'LineWidth', 3);
    xlabel('x');
    ylabel('h(x,t)');
    title(sprintf('Interface Evolution | Time = %.3f', current_time));
    grid on;
    ylim([min_h - margin, max_h + margin]);
    
    % Add progress in title
    sgtitle(sprintf('Frame %d/%d', i, length(sorted_files)));
    
    drawnow;
    pause(0.5); % Adjust speed here
end