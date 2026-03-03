clearvars; close all; clc;

%% 1. Data Acquisition & Sorting
files = dir('newrescaled_h_at_t*.dat');
time_points = zeros(length(files), 1);
for i = 1:length(files)
    t_idx = strfind(files(i).name, 't');
    time_points(i) = str2double(extractBetween(files(i).name, t_idx(2)+1, '.dat'));
end
[sorted_times, sort_idx] = sort(time_points);
sorted_files = files(sort_idx);

%% 2. Peak-Tracking Snapshot Selection
n_snapshots = length(sorted_files); % Higher density for a smoother surface
selected_indices = zeros(n_snapshots, 1);
selected_indices(1) = find(sorted_times > 0, 1); 

current_idx = selected_indices(1);
for i = 2:n_snapshots
    data_curr = readmatrix(sorted_files(current_idx).name);
    [~, max_p] = max(data_curr(:,2));
    z_curr = data_curr(max_p, 1);
    
    found = false;
    search_start = current_idx + 2; % Smaller skip for surface continuity
    for j = search_start:length(sorted_files)
        data_next = readmatrix(sorted_files(j).name);
        [~, p_next] = max(data_next(:,2));
        if abs(data_next(p_next, 1) - z_curr) < 0.5
            selected_indices(i) = j;
            current_idx = j;
            found = true;
            break;
        end
    end
    if ~found, selected_indices(i:end) = []; break; end
end
selected_indices(selected_indices == 0) = [];

%% 3. Surface Reconstruction
% Pre-read grid size
sample = readmatrix(sorted_files(selected_indices(1)).name);
z_vec = sample(:, 1);
Nz = length(z_vec);
Nt = length(selected_indices);

% Build 2D Matrices for surf()
[Z_mesh, T_mesh] = meshgrid(z_vec, sorted_times(selected_indices));
R_mesh = zeros(size(Z_mesh));

for i = 1:Nt
    data = readmatrix(sorted_files(selected_indices(i)).name);
    R_mesh(i, :) = data(:, 2);
end

%% 4. Rendering the Smooth Surface
figure('Position', [100, 100, 1100, 800], 'Color', 'w');

% 'EdgeColor', 'none' makes the surface smooth instead of a mesh
s = surf(Z_mesh, T_mesh, R_mesh, 'EdgeColor', 'none', 'FaceAlpha', 0.9);

% Professional Shading and Lighting
shading interp;      % Interpolates colors across faces for smoothness
lighting gouraud;    % Smooths out the appearance of polygons
camlight right;      % Adds a light source to highlight the wave ridges

% Aesthetics
cb = colorbar;
colormap(jet(256));
grid on; view(-30, 35);
xlabel('Spatial Domain $z$', 'Interpreter', 'latex');
ylabel('Time $t$', 'Interpreter', 'latex');
zlabel('Interface Height $r$', 'Interpreter', 'latex');
title('Smooth Tracked Interface Evolution', 'Interpreter', 'latex');
axis tight;