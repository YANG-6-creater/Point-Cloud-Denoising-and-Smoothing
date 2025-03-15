clc;
clear;
close all;

%% When the value of k remains constant and lamda = 0, 0.5, 1, 1.5, 2, 2.5, 3, the denoising effect of the improved statistical filtering
[FileName, PathName] = uigetfile({'*.asc'; '*.txt'; '*.xlsx'}, 'Select point cloud data');

tic;
str = [PathName FileName];
pointdatas = importdata(str); % Load the point cloud data

%% Add randomly distributed noise with a uniform distribution
a = min(pointdatas, [], 1);     % Minimum values of xyz coordinates
b = max(pointdatas, [], 1);      % Maximum values of xyz coordinates
noise = a + (b - a) .* rand(2000, 3);  % 2000 random noise points, a noise point cloud with coordinates in the range of (a, b)
mu = 0;         % Mean
sigma = 0.11;    % Standard deviation
noises = a + (b - a) .* normrnd(mu, sigma, [500, 3]);  % 500 Gaussian noise points
noises1 = a + (b - a) .* normrnd(mu, 0.08, [500, 3]);  % 500 Gaussian noise points
noises2 = a + (b - a) .* normrnd(mu, 0.02, [500, 3]);  % 500 Gaussian noise points
noises3 = a + (b - a) .* normrnd(mu, 0.05, [500, 3]);  % 500 Gaussian noise points
pointdata = [pointdatas; noise; noises; noises1; noises2; noises3];
point = pointdata';
ptCloud = pointCloud(pointdata);
figure;
pcshow(ptCloud); % Display the point cloud data


% % Improved statistical filtering
k = 5; 
lamda = [0 0.5 1 1.5 2 2.5 3];
[pointdata_rows, pointdata_list] = size(pointdata);
s = [];
for i = 1:7
    final_faces = []; 
    [dist_m, ~, ~, ~] = fajuli(pointdata', k);
    sss = mean(dist_m) + lamda(i) * std(dist_m, 0, 2);
    for ii = 1:pointdata_rows
        point = pointdata(ii, :); 
        if dist_m(ii) < sss
            final_faces = [final_faces; point];
        end   
    end
    s = [s; size(final_faces, 1)];
end

k = 15; 
s1 = [];
for i = 1:7
    final_faces = []; 
    [dist_m, ~, ~, ~] = fajuli(pointdata', k);
    sss = mean(dist_m) + lamda(i) * std(dist_m, 0, 2);
    for ii = 1:pointdata_rows
        point = pointdata(ii, :); 
        if dist_m(ii) < sss
            final_faces = [final_faces; point];
        end   
    end
    s1 = [s1; size(final_faces, 1)];
end

k = 25; 
s2 = [];
for i = 1:7
    final_faces = []; 
    [dist_m, ~, ~, ~] = fajuli(pointdata', k);
    sss = mean(dist_m) + lamda(i) * std(dist_m, 0, 2);
    for ii = 1:pointdata_rows
        point = pointdata(ii, :); 
        if dist_m(ii) < sss
            final_faces = [final_faces; point];
        end   
    end
    s2 = [s2; size(final_faces, 1)];
end

sdata = pointdata;
%% Statistical filtering for denoising
[pointdata_rows, pointdata_list] = size(sdata);

k = 5; 
st = [];
for i = 1:7
    final_d = [];
    m = [];
    for ii = 1:pointdata_rows % Arithmetic sequence, first term 1, last term pointdata_rows, common difference 1.
        point = sdata(ii, :);  % The ii-th row of the matrix.
        [indices, dists] = findNearestNeighbors(ptCloud, point, k);   
        % m = mean(dists); % Calculate the mean
        m = [m; mean(dists)];
    end
    ss2 = std(m, 0, 1); % Calculate the standard deviation    
    dt = mean(m) + lamda(i) * ss2;  
    for ii = 1:pointdata_rows 
        point = sdata(ii, :);  % The ii-th row of the matrix.
        if m(ii) < dt  
            final_d = [final_d; point];
        end
    end
    st = [st; size(final_d, 1)];
end

k = 15;
st1 = [];
for i = 1:7
    final_d = [];
    m = [];
    for ii = 1:pointdata_rows % Arithmetic sequence, first term 1, last term pointdata_rows, common difference 1.
        point = sdata(ii, :);  % The ii-th row of the matrix.
        [indices, dists] = findNearestNeighbors(ptCloud, point, k);  
        m = [m; mean(dists)];
    end
    ss2 = std(m, 0, 1); % Calculate the standard deviation    
    dt = mean(m) + lamda(i) * ss2;  
    for ii = 1:pointdata_rows 
        point = sdata(ii, :);  % The ii-th row of the matrix.
        if m(ii) < dt  
            final_d = [final_d; point];
        end
    end
    st1 = [st1; size(final_d, 1)];
end

k = 25; 
st2 = [];
for i = 1:7
    final_d = [];
    m = [];
    for ii = 1:pointdata_rows % Arithmetic sequence, first term 1, last term pointdata_rows, common difference 1.
        point = sdata(ii, :);  % The ii-th row of the matrix.
        [indices, dists] = findNearestNeighbors(ptCloud, point, k);   
        % m = mean(dists); % Calculate the mean
        m = [m; mean(dists)];
    end
    ss2 = std(m, 0, 1); % Calculate the standard deviation    
    dt = mean(m) + lamda(i) * ss2;  
    for ii = 1:pointdata_rows 
        point = sdata(ii, :);  % The ii-th row of the matrix.
        if m(ii) < dt  
            final_d = [final_d; point];
        end
    end
    st2 = [st2; size(final_d, 1)];
end

figure;
% Resize the figure
set(gcf, 'Position', [100 100 260 220]);
set(gca, 'Position', [.13 .17 .80 .74]);
figure_FontSize = 8;
set(get(gca, 'XLabel'), 'FontSize', figure_FontSize, 'Vertical', 'top');
set(get(gca, 'YLabel'), 'FontSize', figure_FontSize, 'Vertical', 'middle');
set(findobj('FontSize', 10), 'FontSize', figure_FontSize);
set(findobj(get(gca, 'Children'), 'LineWidth', 0.5), 'LineWidth', 2);

y = 0:0.5:3;
p12 = plot(y, st, '--r'); 
hold on
p11 = plot(y, s, '-b');   % z is the average distance difference
hold on
p21 = plot(y, s, '-sk');   % z is the average distance difference
hold on
p22 = plot(y, s1, '-dk');   % z is the average distance difference
hold on
p23 = plot(y, s2, '-*k');   % z is the average distance difference
hold on
p1 = plot(y, s, '-bs');   % z is the average distance difference
hold on
p2 = plot(y, s1, '-bd');   % z is the average distance difference
hold on
p3 = plot(y, s2, '-b*');   % z is the average distance difference
hold on
p4 = plot(y, st, '--rs'); 
hold on
p5 = plot(y, st1, '--rd'); 
hold on
p6 = plot(y, st2, '--r*'); 
hold on
legend('Statistical denoising algorithm', 'Improved statistical denoising algorithm', 'K = 5', 'K = 15', 'K = 25');
hold on;
legend();
xlabel('Standard deviation multiple / \lambda');
ylabel('Post - denoising point cloud count');
hold on;