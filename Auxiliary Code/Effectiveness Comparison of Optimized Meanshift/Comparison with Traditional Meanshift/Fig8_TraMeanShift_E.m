% 1. Select and load point cloud data
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
% Generate Gaussian noise with a mean of 0 and a standard deviation of sigma, with a quantity of 20% of the point cloud
s_noise = normrnd(0, sigma, round(0.2 * size(a, 1)), 3);
% Select a part of the original point cloud every 5 points and add Gaussian noise
s_simulated = a(1:5:5 * round(0.2 * size(a, 1)), 1:3) + s_noise;
% Combine the noisy points and the original point cloud data
ss = [s_simulated; a(:, 1:3)];
% Create a point cloud object ptCloudA containing simulated noisy points
ptCloudA = pointCloud(s_simulated);
% Create a point cloud object ptCloudB containing the original point cloud and simulated noisy points
ptCloudB = pointCloud(ss);

% Display views of the two point cloud objects ptCloud and ptCloudB respectively
figure;
% Set the position and size of the new figure window
set(gcf, 'Position', [100 100 260 220]);
% Set the position and size of the new axes
set(gca, 'Position', [.13 .17 .80 .74]);
% Display the point cloud after mean shift
pcshow(ptCloud);
axis off;

figure;
% Set the position and size of the new figure window
set(gcf, 'Position', [100 100 260 220]);
% Set the position and size of the new axes
set(gca, 'Position', [.13 .17 .80 .74]);
% Display the point cloud after mean shift
pcshow(ptCloudB);
% Turn off the axis display
axis off;

% 4. Set K-nearest neighbor search parameters
% Transpose the combined point cloud data for subsequent K-nearest neighbor search
P = ss';
% Set the number of nearest neighbors for K-nearest neighbor search to 25
k = 25;

% 5. Find the K-nearest neighbors of each point
% Use the knnsearch function to perform K-nearest neighbor search, find the indices of the k nearest neighbors of each point,
% and transpose the result and store it in neighbors
neighbors = transpose(knnsearch(transpose(P), transpose(P), 'k', k + 1));

% 6. Initialize mean shift parameters
% Use the Funpca function to perform principal component analysis on the data P to estimate the normal vector of each point
% Ignore the second output parameter and store the normal vector in pn
[pn, ~] = Funpca(P, k, neighbors);
% Transpose the normal vector matrix
nf = pn';
% Use the combined point cloud data as the point data for subsequent calculations
dataPts = ss;
% Set the bandwidth parameter of the Gaussian kernel function
sigmac = 1.3;
% Set the termination threshold for mean shift iteration
tol = 0.004;
% Get the number of points and dimensions of the point cloud data
[numPts, ~] = size(dataPts);

% 7. Initialize variables
% Initialize a zero matrix z to store the point cloud positions after mean shift
z = zeros(numPts, 3);

% 8. Perform mean shift on each point (still needs modification)
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

    % 9. Calculate the Gaussian kernel
    % Initialize vectors to store the distance differences in direction and position
    distance_f = zeros(1, k);
    % Initialize vectors to store the Gaussian kernel weights based on direction and position differences
    gausskernel_f = zeros(1, k);
    % Initialize vectors to store the position differences
    distance_euclidean = zeros(1, k);
    % Initialize vectors to store the Gaussian kernel weights based on position differences
    gausskernel_euclidean = zeros(1, k);
    for i = 1:k
        % Get the index of the current nearest neighbor point
        b = neighbors(i + 1, n);
        % Calculate the square of the distance difference between the current point and the neighbor point in direction and position
        distance_f(i) = norm((np - nf(b, :))' * (x - dataPts(b, :)) / sigmac).^2;
        % Calculate the Gaussian kernel weight based on direction and position differences
        gausskernel_f(i) = exp(-distance_f(i) / 2);
        % Calculate the square of the position difference between the current point and the neighbor point
        distance_euclidean(i) = norm(x - dataPts(b, :)).^2;
        % Calculate the Gaussian kernel weight based on position differences
        gausskernel_euclidean(i) = exp(-distance_euclidean(i) / 2);
    end

    % 10. Normalize the Gaussian kernel and cluster the neighborhood according to distance differences
    % Sum the Gaussian kernel weights based on position differences
    sumgauss = sum(gausskernel_euclidean);
    % Sum the Gaussian kernel weights based on direction and position differences
    sumgauss_f = sum(gausskernel_f);
    % Normalize the Gaussian kernel weights based on position differences
    ps = gausskernel_euclidean / sumgauss;
    % Normalize the Gaussian kernel weights based on direction and position differences
    ps_f = gausskernel_f / sumgauss_f;
    % Calculate the ratio of the sum of Gaussian kernel weights of direction and position differences to the sum of Gaussian kernel weights of position differences
    mv = sumgauss_f / sumgauss;

    % 11. Calculate mean shift and update the neighborhood
    % Calculate the new position of the current point according to the normalized Gaussian kernel weights based on position differences
    x_d = ps * dataPts(neighbors(2:k + 1, n), :);
    % Calculate the new position of the current point according to the normalized Gaussian kernel weights based on direction and position differences
    x_f = ps_f * dataPts(neighbors(2:k + 1, n), :);

    % Initialize the error to a very large value to ensure entering the iteration loop
    error = realmax;

    % 12. Iteratively update mean shift to smooth the point cloud
    % Continue iterative updating when the error is greater than the set threshold
    while error > tol
        % Save the point position of the previous iteration
        oldx = x_d; % Use the neighborhood points updated based on distance differences
        % Find the indices and distances of the k nearest points to oldx in the point cloud object ptCloudB
        [indices, ~] = findNearestNeighbors(ptCloudB, oldx, k);
        % Select these k nearest neighbor points from the point cloud object ptCloudB
        plo = select(ptCloudB, indices);

        % 13. Calculate the new Gaussian kernel
        % Initialize a vector to store the new distances
        newdistance = zeros(1, k);
        % Initialize a vector to store the new Gaussian kernel weights
        newgausskernel = zeros(1, k);
        % Iterate through these k nearest neighbor points
        for i = 1:k
            % Calculate the square of the distance between the current point and the neighbor point and normalize it
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
        % Update the position of the current point
        x_d = newx;
        % Calculate the error between the new position and the old position
        error = norm(newx - oldx);
    end
    % Store the final point position in the z matrix
    z(n, :) = x_d;
end
ptt1 = pointCloud(z);
figure;
% Set the position and size of the new figure window
set(gcf, 'Position', [100 100 260 220]);
% Set the position and size of the new axes
set(gca, 'Position', [.13 .17 .80 .74]);
% Display the point cloud after mean shift
pcshow(ptt1);
% Turn off the axis display
axis off;