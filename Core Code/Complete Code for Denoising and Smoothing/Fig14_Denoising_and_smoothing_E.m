clc; % Clear the command window
clear; % Clear variables in the workspace
close all; % Close all figure windows
%% Select and load the point cloud data
[FileName, PathName] = uigetfile({'*.asc';'*.txt';'*.xlsx'}, 'Select point cloud data'); % File selection dialog box
str = [PathName FileName]; % Construct the full file path
pointdatas_bef = importdata(str); % Load the point cloud data (supports asc/txt/xlsx formats)
%% Display the point cloud before processing
ptCloud = pointCloud(pointdatas_bef); % Create a point cloud object
figure; % Create a new figure window
pcshow(ptCloud); % Visualize and display the point cloud
axis off; % Hide the coordinate axes
title('Point cloud before denoising and smoothing');
disp(['Total number of initial point cloud points (including noise): ', num2str(size(pointdatas_bef, 1))]); % Display the total number of points
%% Improved statistical filtering denoising algorithm
pointdata = pointdatas_bef; % Copy the noisy point cloud data
k = 30;     % Maximum number of neighboring points (parameter can be adjusted)
lambda = 0.5; % Three lambda thresholds (parameter can be adjusted)
[pointdata_rows, ~] = size(pointdata); % Get the data dimensions
% Perform filtering processing
final_faces_denoising = []; % Store the filtered point cloud
 % Calculate the average distance from each point to its k-nearest neighbors (requires a custom function fajuli)
[dist_m, ~, ~, ~] = fajuli(pointdata', k); % Transpose and pass the parameters
% Dynamic threshold calculation (d_y = mean + lambda*standard deviation)
sss = mean(dist_m) + lambda * std(dist_m, 0, 2); % Calculate the standard deviation horizontally
for ii = 1:pointdata_rows
    point = pointdata(ii, :); % Extract a single point
    if dist_m(ii) < sss % If the distance is less than the threshold, keep it
        final_faces_denoising = [final_faces_denoising; point]; % Vertically concatenate valid points
    end
end
% Visualize the denoising effect of filtering
figure;
ptCloudss_denoising=pointCloud(final_faces_denoising);
pcshow(ptCloudss_denoising);
axis off;
%% Improved Meanshift point cloud smoothing
%% K-nearest neighbor search (preprocessing)
ss=final_faces_denoising;
% Create a point cloud object to be smoothed
ptCloudB = pointCloud(ss);
P = ss';     % Transpose the point cloud matrix


k = 25;      % Set the K value
neighbors = transpose(knnsearch(transpose(P), transpose(P), 'k', k+1)); 
% Structure of the neighbors matrix: The first column is the point index, and the subsequent k columns are the nearest neighbor indices
%% Comprehensive noise degree calculation (call a custom function)
[hunhe_d, qiuzaosheng, bianfen, dist_m, dist, pn, pw] = funzhong_zaoshengdu(P, k, neighbors);
nf = pn';    % Transpose the normal vector matrix
%% Mean shift parameter initialization
% Use principal component analysis to obtain the normal vector (ignore secondary outputs)
[pn, ~] = Funpca(P, k, neighbors);
nf = pn';     % Update the normal vector matrix
dataPts = ss;  % Use the noisy point cloud for subsequent processing
sigmac = 0.3;  % Gaussian kernel bandwidth parameter
tol = 0.004;   % Iteration termination threshold
[numPts, ~] = size(dataPts); % Get the dimensional information of the point cloud
%% Variable initialization
z = zeros(numPts, 3); % Store the mean shift results
%% Core algorithm of mean shift
for n = 1:numPts
    x = dataPts(n, :);     % Coordinates of the current point
    np = nf(n, :);         % Normal vector of the current point
    
    % Get the normal vector matrix of the k-nearest neighbors
    npm = nf(neighbors(:,n), :);
    
    % Use the improved mean shift algorithm (including the regularization term)
    [clustCent, data2cluster, ~] = meanshift(npm', 0.4);
    fa = clustCent(:, data2cluster(1)); % Select the regularized vector of the cluster center
      %% Gaussian kernel calculation (combined weight of direction and position)
    distance_f = zeros(1,k);       % Square of the combined distance of direction and position
    gausskernel_f = zeros(1,k);     % Combined kernel weight
    distance_euclidean = zeros(1,k);% Square of the Euclidean distance
    gausskernel_euclidean = zeros(1,k);% Kernel weight of position alone
    
    for i = 1:k
        b = neighbors(i+1, n);     % Index of the i-th neighbor
        % Calculate the combined distance of direction and position (after normalization)
        distance_f(i) = norm( (np - nf(b,:))' * (x - dataPts(b,:)) / sigmac )^2;
        gausskernel_f(i) = exp(-distance_f(i)/2);
        
        % Calculate the distance of position alone (after normalization)
        distance_euclidean(i) = norm( (x - dataPts(b,:)) / sigmac )^2;
        gausskernel_euclidean(i) = exp(-distance_euclidean(i)/2);
    end
    %% Normalize the kernel function
    sumgauss = sum(gausskernel_euclidean); % Sum of the position kernel
    sumgauss_f = sum(gausskernel_f);     % Sum of the combined kernel
    ps = gausskernel_euclidean / sumgauss;   % Normalized position kernel
    ps_f = gausskernel_f / sumgauss_f;     % Normalized combined kernel
    mv = sumgauss_f / sumgauss;           % Kernel weight scaling factor
       %% Mean shift estimation
    x_d = ps * dataPts(neighbors(2:k+1,n),:); % Weighted average of the position kernel
    x_f = ps_f * dataPts(neighbors(2:k+1,n),:); % Weighted average of the combined kernel
    
    error = realmax; % Initialize the maximum error
    %% Iterative optimization
    while error > tol
        oldx = x_f; % Result of the previous iteration
        
        % Find the k-nearest neighbors in the mixed point cloud
        [indices, ~] = findNearestNeighbors(ptCloudB, oldx, k);
        plo = select(ptCloudB, indices); % Get the nearest neighbor point cloud data
        %% Calculate the new kernel function (based on the real point cloud)
        newdistance = zeros(1,k);
        newgausskernel = zeros(1,k);
        for i = 1:k
            newdistance(i) = norm( (oldx - plo.Location(i,:)) / sigmac )^2;
            newgausskernel(i) = exp(-newdistance(i)/2);
        end
        
        sumnewgauss = sum(newgausskernel);
        newp = newgausskernel / sumnewgauss;
        x_m = newp * plo.Location; % Calculate the new position
        
        % Update formula (including the noise suppression term)
        x_f = x + hunhe_d(n) * dot(np, (x_m - x)) * np;
        error = norm(x_m - oldx); % Calculate the error
        end
    zz(n,:) = x_f; % Save the optimized coordinates
end
%% Visualize the results
ptt2 = pointCloud(zz); % Construct the denoised point cloud
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
pcshow(ptt2);
axis off;