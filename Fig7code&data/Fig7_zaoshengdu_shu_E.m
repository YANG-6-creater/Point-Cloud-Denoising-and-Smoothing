% Clear the command window, workspace, and close all figures
clc;
clear;
close all;

% Prompt the user to select point cloud data
[FileName, PathName] = uigetfile({'*.asc'; '*.txt'; '*.xlsx'}, 'Select point cloud data');

% Start the timer
tic;

% Combine the path and file name
str = [PathName FileName];

% Load the point cloud data
pointdatas = importdata(str); 

%% Add Gaussian noise
ptCloud = pointCloud(pointdatas);
a = ptCloud.Location;
sigma = 0.008;
% Add noise with mean parameter mu = 0 and standard deviation parameter sigma = 0.008
s_noise = normrnd(0, sigma, round(0.2 * size(a, 1)), 3);
s_simulated = a(1:5:5 * round(0.2 * size(a, 1)), 1:3) + s_noise;

% pointdata = s_simulated;
ss = [s_simulated; a(:, 1:3)];
pointdata = ss;
ptCloudA = pointCloud(s_simulated); % Simulated Gaussian noise points near the outer surface of the rabbit
ptCloudB = pointCloud(ss); % Simulated Gaussian noise points near the outer surface and the original rabbit point cloud

% Set the figure properties
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
set(gca, 'tickdir', 'in');
figure_FontSize = 8;

% Display the simulated noise points and the original point cloud in one window
pcshowpair(ptCloudA, ptCloud); 
axis off;

% Set the number of nearest neighbors
k = 25;

%%
[qiuzaosheng, bianfen, dist_m, dist, n, w] = fun_zaoshengdu(pointdata', k);
[pointdata_rows, pointdata_list] = size(pointdata);
ping_d = dist ./ (dist + dist_m); % Planar fitting noise degree d/(d_m + d)
ping_d = ping_d'; % Calculate the noise probability of each point under planar fitting and store it as a feature vector

final_face = []; % Minimum ball noise degree
for ii = 1:pointdata_rows
    point = pointdata(ii, :);
    [indices, dists] = findNearestNeighbors(ptCloudB, point, k);
    p = select(ptCloudB, indices);
    p = p.Location;
    center = mean(p, 1); % Calculate the center coordinates of the current point
    center_dist = pdist2(point, center); % Calculate the Euclidean distance from the current point to the center point
    center_dists = pdist2(point, p); % Calculate the Euclidean distances from the current point to its nearest neighbors
    qiu_d = center_dist / (2 * mean(center_dists) / sqrt(k) + center_dist); % Reflect the degree of the current point deviating from the neighborhood center
    final_face = [final_face; qiu_d]; % Calculate the noise probability of each point under the minimum bounding ball and store it as a feature vector
end

hunhe_d = 0.4 * ping_d + 0.6 * final_face;
pingyuzhi = 0.68; % Planar fitting noise probability threshold
qiuyuzhi = 0.38; % Minimum bounding ball fitting noise probability threshold
Hyuzhi = 0.50;    % Mixed fitting noise probability threshold

ping_p = [];
qiu_p = [];
h_p = [];
for iis = 1:pointdata_rows
    point = pointdata(iis, :);
    if ping_d(iis) > pingyuzhi
        ping_p = [ping_p; point];
    end
    if final_face(iis) > qiuyuzhi
        qiu_p = [qiu_p; point];
    end
    if hunhe_d(iis) > Hyuzhi
        h_p = [h_p; point];
    end
end

% Create a new figure and display the Gaussian noise identified by planar fitting probability and the original point cloud
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
set(gca, 'tickdir', 'in');
figure_FontSize = 8;
ptCloudss = pointCloud(ping_p);
pcshowpair(ptCloudss, ptCloud); 
axis off;

% Create a new figure. Note: There was a variable 'ptClouds' which might be a typo.
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
set(gca, 'tickdir', 'in');
figure_FontSize = 8;
% Assume you meant to use 'qiu_p' here to create the point cloud
ptClouds = pointCloud(qiu_p); 
pcshowpair(ptClouds, ptCloud); % Display the Gaussian noise identified by minimum bounding ball noise probability and the original point cloud
axis off;

% Create a new figure and display the Gaussian noise identified by mixed probability and the original point cloud
figure;
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
set(gca, 'tickdir', 'in');
figure_FontSize = 8;
ptCloudsh = pointCloud(h_p);
pcshowpair(ptCloudsh, ptCloud);
axis off;