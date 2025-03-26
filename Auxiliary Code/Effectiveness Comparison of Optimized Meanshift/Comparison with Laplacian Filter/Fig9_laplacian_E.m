% Point cloud smoothing algorithm: Laplacian smoothing
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

% 4. Implement the Laplacian smoothing algorithm
fprintf('Starting Laplacian smoothing...\n');

% Get the point cloud data
points = ptCloudB.Location;

% Laplacian smoothing parameters
k = 25;           % Number of neighboring points
iterations = 3;  % Number of iterations
lambda = 0.5;     % Smoothing coefficient

% Build a KD tree to accelerate neighborhood search
kdtree = KDTreeSearcher(points);

% Initialize the result
laplacian_smoothed = points;

% Iterative smoothing process
for iter = 1:iterations
    % Create new point cloud coordinates
    new_points = zeros(size(laplacian_smoothed));
    
    % Process each point
    parfor i = 1:size(laplacian_smoothed, 1)
        % Find the k nearest neighbors of each point
        [indices, ~] = knnsearch(kdtree, laplacian_smoothed(i,:), 'K', k+1);
        indices = indices(2:end);  % Exclude the point itself
        
        % Calculate the centroid of the neighborhood
        neighbor_points = laplacian_smoothed(indices, :);
        centroid = mean(neighbor_points, 1);
        
        % Calculate the Laplacian coordinates
        laplacian_vector = centroid - laplacian_smoothed(i,:);
        
        % Apply smoothing
        new_points(i,:) = laplacian_smoothed(i,:) + lambda * laplacian_vector;
    end
    
    % Update the point cloud coordinates
    laplacian_smoothed = new_points;
    
    % Display the progress
    fprintf('Laplacian smoothing: Completed iteration %d/%d\n', iter, iterations);
end

% Create a point cloud object after Laplacian smoothing
ptCloudLaplacian = pointCloud(laplacian_smoothed);

% Display the result of Laplacian smoothing
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
pcshow(ptCloudLaplacian);
title('Result of Laplacian Smoothing');
axis off;

% Save the result of Laplacian smoothing
pcwrite(ptCloudLaplacian, 'laplacian_smoothed.ply');
fprintf('Laplacian smoothing completed. The result has been saved as laplacian_smoothed.ply\n');