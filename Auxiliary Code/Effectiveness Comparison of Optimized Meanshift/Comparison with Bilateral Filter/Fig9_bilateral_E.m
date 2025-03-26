% Point cloud smoothing algorithm: Bilateral filtering
% Implemented based on Matlab 2024a
% Process the point cloud data in ptCloudB

% 1. Select and load point cloud data
clc;
clear;
close all;
[FileName, PathName] = uigetfile({'*.asc'; '*.txt'; '*.xlsx'}, 'Select point cloud data');

str = [PathName FileName];
pointdatas = importdata(str);
% 2. Create a point cloud object
ptCloud = pointCloud(pointdatas);

% 3. Display the point cloud data
a = ptCloud.Location;
n = 0.001;
s_noise = normrnd(0, n, round(0.2 * size(a, 1)), 3);
s_simulated = a(1:5:5 * round(0.2 * size(a, 1)), 1:3) + s_noise;
ss = [s_simulated; a(:, 1:3)];
ptCloudA = pointCloud(s_simulated); 
ptCloudB = pointCloud(ss); 

% Display the original point cloud
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
pcshow(ptCloud);
title('Original Point Cloud');
axis off;

% Display the noisy point cloud
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
pcshow(ptCloudB);
title('Noisy Point Cloud');
axis off;

% 5. Implement the bilateral filtering algorithm
fprintf('Starting bilateral filtering...\n');

% Bilateral filtering parameters
k_bilateral = 25;         % Number of neighboring points
sigma_s = 0.05;           % Spatial domain standard deviation
sigma_r = 0.05;           % Range domain standard deviation
iterations_bilateral = 3; % Number of iterations

% Get the point cloud data
points = ptCloudB.Location;

% Calculate the point cloud normals (for bilateral filtering)
normals = pcnormals(ptCloudB);

% Build a KD tree
kdtree = KDTreeSearcher(points);

% Initialize the result
bilateral_smoothed = points;

% Iterative smoothing process
for iter = 1:iterations_bilateral
    % Create new point cloud coordinates
    new_points = zeros(size(bilateral_smoothed));
    
    % Process each point
    parfor i = 1:size(bilateral_smoothed, 1)
        % Find the k nearest neighbors of each point
        [indices, distances] = knnsearch(kdtree, bilateral_smoothed(i,:), 'K', k_bilateral+1);
        indices = indices(2:end);  % Exclude the point itself
        distances = distances(2:end);
        
        % Get the neighboring points
        neighbor_points = bilateral_smoothed(indices, :);
        neighbor_normals = normals(indices, :);
        current_normal = normals(i, :);
        
        % Calculate the spatial domain weights (based on distance)
        weights_s = exp(-distances.^2 / (2 * sigma_s^2));
        
        % Calculate the range domain weights (based on normal vector differences)
        proj_dists = zeros(length(indices), 1);
        for j = 1:length(indices)
            % Calculate the distance from the point to the tangent plane (projection along the normal vector)
            vector_diff = neighbor_points(j,:) - bilateral_smoothed(i,:);
            proj_dists(j) = abs(dot(vector_diff, current_normal));
        end
        weights_r = exp(-proj_dists.^2 / (2 * sigma_r^2));
        
        % Combine the two weights
        weights = weights_s .* weights_r;
        weights = weights / sum(weights);
        
        % Calculate the weighted average position
        if sum(weights) > 0
            weighted_shift = zeros(1, 3);
            for j = 1:length(indices)
                weighted_shift = weighted_shift + weights(j) * (neighbor_points(j,:) - bilateral_smoothed(i,:));
            end
            new_points(i,:) = bilateral_smoothed(i,:) + weighted_shift;
        else
            new_points(i,:) = bilateral_smoothed(i,:);
        end
    end
    
    % Update the point cloud coordinates
    bilateral_smoothed = new_points;
    
    % Display the progress
    fprintf('Bilateral filtering: Completed iteration %d/%d\n', iter, iterations_bilateral);
end

% Create a point cloud object after bilateral filtering
ptCloudBilateral = pointCloud(bilateral_smoothed);

% Display the result of bilateral filtering
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
pcshow(ptCloudBilateral);
title('Result of Bilateral Filtering');
axis off;

% Save the result of bilateral filtering
pcwrite(ptCloudBilateral, 'bilateral_smoothed.ply');
fprintf('Bilateral filtering completed. The result has been saved as bilateral_smoothed.ply\n');