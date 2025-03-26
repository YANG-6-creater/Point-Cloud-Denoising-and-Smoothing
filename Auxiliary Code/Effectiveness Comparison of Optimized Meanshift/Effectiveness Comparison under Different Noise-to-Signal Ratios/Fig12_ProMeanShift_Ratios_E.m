%1. Select and load point cloud data
clc;
clear;
close all;

% Select a point cloud data file through a file dialog, supporting .asc, .txt, and .xlsx formats
[FileName, PathName] = uigetfile({'*.asc'; '*.txt'; '*.xlsx'}, 'Select point cloud data');
str = [PathName FileName];
pointdatas = importdata(str);

% 2. Create a point cloud object
% Create a point cloud object ptCloud using the imported point cloud data
ptCloud = pointCloud(pointdatas);

% 3. Display the point cloud data
% Extract the position data from the point cloud object
a = ptCloud.Location;
% Set the standard deviation of Gaussian noise
sigma = 0.001;
% Generate Gaussian noise with a mean of 0, a standard deviation of sigma, and a quantity of 20% of the point cloud
s_noise = normrnd(0, sigma, round(0.20 * size(a, 1)), 3);
% Select a part of the original point cloud every 5 points and add Gaussian noise
s_simulated = a(1:5:5 * round(0.20 * size(a, 1)), 1:3) + s_noise;
% Combine the noisy points and the original point cloud data
ss = [s_simulated; a(:, 1:3)];
% Create a point cloud object ptCloudA containing simulated noisy points
ptCloudA = pointCloud(s_simulated);
% Create a point cloud object ptCloudB containing the original point cloud and simulated noisy points
ptCloudB = pointCloud(ss);

% Display views of the two point cloud objects ptCloudA and ptCloudB respectively
figure;
% Set the position and size of the new figure window
set(gcf, 'Position', [100 100 260 220]);
% Set the position and size of the new axes
set(gca, 'Position', [.13 .17 .80 .74]);
% Display the original point cloud
pcshow(ptCloud);
title('Original Point Cloud');
pcwrite(ptCloud, 'ptCloud.ply');
axis off;

figure;
% Set the position and size of the new figure window
set(gcf, 'Position', [100 100 260 220]);
% Set the position and size of the new axes
set(gca, 'Position', [.13 .17 .80 .74]);
% Display the noisy point cloud
pcshow(ptCloudB);
title('Noisy Point Cloud');
% Turn off the axis display
axis off;

P = ss';
k = 25;
% Perform K-nearest neighbor search and transpose the result.
% neighbors is an index matrix. The first row represents the point number, and the next k rows represent the K-nearest neighbor points.
% It records each point and its surrounding k points.
neighbors = transpose(knnsearch(transpose(P), transpose(P), 'k', k + 1));

% Calculate comprehensive noise degree and related parameters
[hunhe_d, qiuzaosheng, bianfen, dist_m, dist, pn, pw] = funzhong_zaoshengdu(P, k, neighbors);
nf = pn';

% 6. Find the K-nearest neighbors of each point
% Use the knnsearch function to perform K-nearest neighbor search, find the indices of the k nearest neighbors of each point,
% and transpose the result and store it in neighbors
neighbors = transpose(knnsearch(transpose(P), transpose(P), 'k', k + 1));

% 7. Initialize mean shift parameters
% Use the Funpca function to perform principal component analysis on the data P to estimate the normal vector of each point
% Ignore the second output parameter and store the normal vector in pn
[pn, ~] = Funpca(P, k, neighbors);
% Transpose the normal vector matrix
nf = pn';
% Use the combined point cloud data as the point data for subsequent calculations
dataPts = ss;
% Set the bandwidth parameter of the Gaussian kernel function
sigmac = 0.3;
% Set the termination threshold for mean shift iteration
tol = 0.004;
% Get the number of points and dimensions of the point cloud data
[numPts, ~] = size(dataPts);

% 8. Initialize variables
% Initialize a zero matrix z to store the point cloud positions after mean shift
z = zeros(numPts, 3);

% 9. Perform mean shift on each point (still needs modification)
% Iterate through each point in the point cloud data
for n = 1:numPts
    % Get the coordinates of the current point
    x = dataPts(n, :);
    % Get the normal vector of the current point
    np = nf(n, :);
    % Get the normal vectors of the current point and its k nearest neighbors
    npm = nf(neighbors(:, n), :);
    % Set the bandwidth for mean shift clustering
    bandWidth = 0.4;
    % Use the meanshift function to cluster the normal vectors of the current point and its neighbors
    % Get the cluster centers, the cluster indices of each point, and the number of clusters
    [clustCent, data2cluster, ~] = meanshift(npm', bandWidth);
    % Replace the original normal vector with the center normal vector of the cluster to which the current point belongs to enhance noise resistance
    fa = clustCent(:, data2cluster(1));

    % 10. Calculate the Gaussian kernel
    % Initialize vectors to store the distance differences in direction and position
    distance_f = zeros(1, k);
    % Initialize vectors to store the Gaussian kernel weights based on direction and position differences
    gausskernel_f = zeros(1, k);

    for i = 1:k
        % Get the index of the current nearest neighbor point
        b = neighbors(i + 1, n);
        % Calculate the square of the difference in direction and position between the current point and the neighbor point
        distance_f(i) = norm((np - nf(b, :))' * (x - dataPts(b, :)) / sigmac).^2;
        % Calculate the Gaussian kernel weight based on direction and position differences
        gausskernel_f(i) = exp(-distance_f(i) / 2);
    end

    % 11. Normalize the Gaussian kernel
    % Sum the Gaussian kernel weights based on direction and position differences
    sumgauss_f = sum(gausskernel_f);
    % Normalize the Gaussian kernel weights based on direction and position differences
    ps_f = gausskernel_f / sumgauss_f;

    % 12. Calculate mean shift and update the neighborhood
    % Calculate the new position of the current point according to the normalized Gaussian kernel weights based on direction and position differences
    x_f = ps_f * dataPts(neighbors(2:k + 1, n), :);

    % Initialize the error to a very large value to ensure entering the iteration loop
    error = realmax;

    % 13. Iteratively update mean shift to smooth the point cloud
    % Continue iterative updating when the error is greater than the set threshold
    while error > tol
        % Save the point position of the previous iteration
        oldx = x_f; % Use the neighborhood points updated based on direction and position differences
        % Find the indices and distances of the k nearest points to oldx in the point cloud object ptCloudB
        [indices, ~] = findNearestNeighbors(ptCloudB, oldx, k);
        % Select these k nearest neighbor points from the point cloud object ptCloudB
        plo = select(ptCloudB, indices);

        % 14. Calculate the new Gaussian kernel
        % Initialize a vector to store the new distances
        newdistance = zeros(1, k);
        % Initialize a vector to store the new Gaussian kernel weights
        newgausskernel = zeros(1, k);
        % Iterate through these k nearest neighbor points
        for i = 1:k
            % Calculate the square of the distance between the clustered point and the neighbor point and normalize it
            newdistance(i) = norm((oldx - plo.Location(i, :)) / sigmac).^2;
            % Calculate the Gaussian kernel weight based on the new distance
            newgausskernel(i) = exp(-newdistance(i) / 2);
        end

        % 14. Normalize the new Gaussian kernel
        % Sum the new Gaussian kernel weights
        sumnewgauss = sum(newgausskernel);
        % Normalize the new Gaussian kernel weights
        newp = newgausskernel / sumnewgauss;

        % 15. Update the position
        % Calculate the new position of the current point according to the normalized new Gaussian kernel weights
        newx = newp * plo.Location;
        x_m = newx;
        % Update the position of the current point
        x_f = x + hunhe_d(n).*dot(np,(x_m - x))*np; % New point
        % Calculate the error between the new position and the old position
        error = norm(newx - oldx);
    end
    % Store the final point position in the zz matrix
    zz(n, :) = x_f;
end

% 16. Display the point cloud after mean shift
ptt2 = pointCloud(zz);

figure;
% Set the position and size of the new figure window
set(gcf, 'Position', [100 100 260 220]);
% Set the position and size of the new axes
set(gca, 'Position', [.13 .17 .80 .74]);
% Display the point cloud after mean shift
pcshow(ptt2);
title('Improved Meanshift Smoothing Result');
% Turn off the axis display
axis off;
% Save the result of the improved MeanShift smoothing
pcwrite(ptt2, 'ProMeanShift_smoothed0.20.ply');
fprintf('Improved MeanShift smoothing completed. The result has been saved as ProMeanShift_smoothed0.20.ply\n');