clc;
clear;

files = dir('rescaled_h_at_t*.dat');
file2 = dir('rescaled_h_at_t*.dat');
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

% Extract and sort the second set of files by time
time_points2 = zeros(1, length(file2));
for i = 1:length(file2)
    filename = file2(i).name;
    % Extract time from filename 
    t_indices = strfind(filename, 't');
    time_str = extractBetween(filename, t_indices(2)+1, '.dat');
    time_points2(i) = str2double(time_str);
end
[sorted_times2, sort_idx2] = sort(time_points2);
sorted_files2 = file2(sort_idx2);

% Initialize times array
times = zeros(1, min(length(sorted_files1), length(sorted_files2)));

% Use the shorter length to avoid index issues
num_files = min(length(sorted_files1), length(sorted_files2));

for n = num_files:-plot_n:1
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
    times(n) = sorted_times1(n);

    % Construct full filename for second dataset
    filename2 = sorted_files2(n).name;

    % Read the data (two-column format: x y)
    data2 = readmatrix(filename2);

    % Check if data is valid and has 2 columns
    if isempty(data2) || size(data2, 2) ~= 2
        fprintf('Warning: %s is empty or has incorrect format. Skipping.\n', filename2);
        continue;
    end

    % Extract coordinates
    x2 = data2(:,1);  % Scale axial coordinate to [0, 2π] for MATLAB comparison
    y2 = data2(:,2);         % Radial coordinate (r)

    
    figure;%('Name', sprintf('Interface at t = %.2f', t));
    plot(x1,y1, 'b.', 'MarkerSize', 10);
    hold on
    plot(x2,y2, 'r-');
    hold off
    xlabel('z (scaled to [0, 2\pi])');
    ylabel('r');
    title(sprintf('Interface at t = %.2f', times(n)));
    grid on;
end