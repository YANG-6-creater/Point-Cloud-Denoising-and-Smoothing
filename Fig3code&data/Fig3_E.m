clc;
clear;
close all;

% Select and load the point cloud data
[FileName, PathName] = uigetfile({'*.asc';'*.txt';'*.xlsx'}, 'Select point cloud data');
tic;
str = [PathName FileName];
pointdatas = importdata(str); % Load the point cloud data

%% Add random noise with uniform and Gaussian distributions
a = min(pointdatas, [], 1); % Minimum values of xyz coordinates
b = max(pointdatas, [], 1); % Maximum values of xyz coordinates

% Generate random noise
noise = a + (b - a) .* rand(2000, 3);          % 2000 uniformly distributed noise points
mu = 0;                                        % Mean
sigma1 = 0.11;                                 % Standard deviation 1
sigma2 = 0.08;                                 % Standard deviation 2
noises = a + (b - a) .* normrnd(mu, sigma1, [500, 3]); % 500 Gaussian noise points
noises1 = a + (b - a) .* normrnd(mu, sigma2, [500, 3]);% 500 Gaussian noise points

% Merge the original point cloud and the noise
sdata = [pointdatas; noise; noises; noises1];

% Display the point cloud with added noise
ptCloud = pointCloud(sdata);
figure;
pcshow(ptCloud);
axis off;
title('Point cloud with random noise');
disp(['Initial number of point cloud points (including noise): ', num2str(size(sdata, 1))]);

%% Improved statistical filtering for denoising
pointdata = sdata;
k = 15;     % Maximum number of neighboring points
ji = 5;     % Step size
lambdas = [0.8, 0.5, 0.05]; % List of different λ values

[pointdata_rows, ~] = size(pointdata);
ss = [];

% Iterate through different λ values for filtering
for lambda = lambdas
    for i = 5:ji:k
        final_faces = [];

        % Calculate the average distance from each point to its k-nearest neighbors (using custom function fajuli)
        [dist_m, ~, ~, ~] = fajuli(pointdata', i);

        % Calculate the distance threshold (d_y = d_ + λ * σ)
        sss = mean(dist_m) + lambda * std(dist_m, 0, 2);

        % Filter and keep points based on the threshold
        for ii = 1:pointdata_rows
            point = pointdata(ii, :);
            if dist_m(ii) < sss
                final_faces = [final_faces; point];
            end
        end

        % Record the number of point cloud points
        ss = [ss; size(final_faces, 1)];

        % Display the point cloud after each filtering
        figure;
        ptCloudss = pointCloud(final_faces);
        pcshow(ptCloudss);
        axis off;
        title(['k = ', num2str(i), ', \lambda = ', num2str(lambda), ' Processed point cloud']);

        % Output the current number of point cloud points
        disp(['Number of point cloud points when k = ', num2str(i), ', λ = ', num2str(lambda), ': ', num2str(size(final_faces, 1))]);
    end
end