If you use this code and dataset, please cite it in the following format:

Article Title: [Enhanced Point Cloud Denoising and Smoothing via Improved Statistical Denoising Algorithm and Global Noise Probability-Based Meanshift Algorithm]

Journal Name: （The Visual Computer）

The core code consists mainly of three parts. Each part will be introduced separately.

The first part:  Fig3_E.

The second part: Fig7_zaoshengdu_shu_E.

The third part:  Fig8_ProMeanShift_E.

Part1 Fig3_E

1.Project Overview

This project provides core code for point cloud denoising, mainly implementing the functions of loading point cloud data, adding random noise (uniform distribution and Gaussian distribution), and using an improved statistical filtering algorithm for denoising. Through different parameter settings, the noisy point cloud is filtered multiple times, and the results after each processing are visualized.

2.Dependencies and Requirements

2.1 Software Environment
 
MATLAB: This code is written and run in the MATLAB environment. It is recommended to use MATLAB R2022
 or higher versions to ensure that the functions in the code can be used properly.
	
2.2 Dependent Functions
 
A custom function fajuli is called in the code, which is used to calculate the average distance from each point to its k - nearest neighbors. Please ensure that the fajuli function file exists in the MATLAB search path or is in the same directory as this code file.

2.3 Data Requirements
 
The code supports loading point cloud data files in .asc, .txt, and .xlsx formats. When running the code, a file selection dialog box will pop up, and you need to select the corresponding point cloud data file.

3.Description and Implementation of Key Algorithms
	
3.1 Adding Random Noise
 
Principle: To simulate the noise in the real scene, the code adds random noise with uniform distribution and Gaussian distribution to the original point cloud data.
Implementation Steps:
Read the original point cloud data and obtain the minimum and maximum values of its xyz coordinates.
Generate 2000 uniformly distributed noise points within the coordinate range of the original point cloud.
Generate 500 Gaussian - distributed noise points with standard deviations of 0.11 and 0.05 respectively, and also limit their range within the coordinate range of the original point cloud.
Merge the original point cloud data and the generated noise points to obtain the noisy point cloud data.
	
3.2 Improved Statistical Filtering for Denoising
 
Principle: Statistical filtering is a commonly used method for point cloud denoising. Its basic idea is to calculate the distance from each point to its neighboring points, and judge whether the point is a noise point according to the statistical information (such as mean and standard deviation) of the distance. Based on the traditional statistical filtering, this code introduces different λ values to adjust the distance threshold to improve the denoising effect.
Implementation Steps:
Initialize parameters, including the maximum number of neighboring points k, the step size ji, and the list of different λ values lambdas.
For each λ value, starting from 5 and increasing to k with the step size ji, perform the filtering process successively.
For each k value, call the fajuli function to calculate the average distance from each point to its k - nearest neighbors.
Calculate the distance threshold d_y = d_ + λ * σ according to the average distance, where d_ is the mean value of the average distance, and σ is the standard deviation.
Traverse each point. If the average distance from the point to its neighboring points is less than the threshold, keep the point; otherwise, regard it as a noise point and remove it.
Record the number of points in the point cloud after each filtering, and visualize the processed point cloud.

4.Code Structure and Running Process


4.1 Code Structure
 
Data Loading: Use the uigetfile function to select and load the point cloud data file, and then use the importdata function to read the data.
Noise Addition: Calculate the minimum and maximum values of the coordinates of the original point cloud, generate uniformly distributed and Gaussian - distributed noise points, and merge them with the original point cloud.
Denoising Processing: Through nested loops, perform multiple filtering processes for different λ values and k values. After each processing, record the number of points and visualize the results.

4.2 Running Process
 
Run the code, and a file selection dialog box will pop up. Select the point cloud data file to be processed.
The code will automatically load the data and add random noise, and display the noisy point cloud.
Start the denoising process. Perform multiple filtering according to different λ values and k values. After each filtering, the processed point cloud will be displayed, and the number of points in the current point cloud will be output in the command window.

Part2 Fig7_zaoshengdu_shu_E

1.Dependencies and Requirements

1.1 Software Environment

MATLAB: The code is written in MATLAB. You need to install the MATLAB software. It is recommended to use a relatively new version (such as R2020b or higher) to ensure the compatibility and performance of the code.

1.2 Toolboxes

Point Cloud Toolbox: The code uses functions like pointCloud, pcshowpair, and findNearestNeighbors, which belong to MATLAB's Point Cloud Toolbox. Before running the code, make sure this toolbox is correctly installed and available.

1.3 Data Requirements

The supported point - cloud data file formats are .asc, .txt, and .xlsx. When running the code, users need to select point - cloud data files in these formats.

2.Key Algorithm Description and Implementation

2.1 Data Loading

Function: Prompt the user to select a point - cloud data file and load it into the MATLAB workspace.
Implementation:
matlab
[FileName, PathName] = uigetfile({'*.asc'; '*.txt'; '*.xlsx'}, 'Select point cloud data');

str = [PathName FileName];

pointdatas = importdata(str); 

2.2 Adding Gaussian Noise

Function: Add Gaussian noise to the point - cloud data to simulate possible noise situations in practical applications.
matlab

ptCloud = pointCloud(pointdatas);

a = ptCloud.Location;

sigma = 0.008;

s_noise = normrnd(0, sigma, round(0.2 * size(a, 1)), 3);

s_simulated = a(1:5:5 * round(0.2 * size(a, 1)), 1:3) + s_noise;

ss = [s_simulated; a(:, 1:3)];

pointdata = ss;

ptCloudA = pointCloud(s_simulated); 

ptCloudB = pointCloud(ss); 

2.3 Noise Level Calculation

3.Noise Level Calculation

3.1 Planar Fitting Noise Level

Function: Evaluate the planar fitting noise level of a point by calculating the distance from the point to its neighborhood plane.
Implementation:
matlab

[qiuzaosheng, bianfen, dist_m, dist, n, w] = fun_zaoshengdu(pointdata', k);

ping_d = dist ./ (dist + dist_m); 

ping_d = ping_d'; 

3.2 Minimum Bounding Sphere Noise Level

Function: Evaluate the minimum bounding sphere noise level of a point by calculating the relationship between the distance from the point to its neighborhood center and the average distance to its neighborhood points.
Implementation:
matlab

final_face = [];

for ii = 1:pointdata_rows
    
		point = pointdata(ii, :);
   
		[indices, dists] = findNearestNeighbors(ptCloudB, point, k);
   
		p = select(ptCloudB, indices);
  
		p = p.Location;
  
		center = mean(p, 1);
  
		center_dist = pdist2(point, center);
   
		center_dists = pdist2(point, p);
   
		qiu_d = center_dist / (2 * mean(center_dists) / sqrt(k) + center_dist);
  
		
	  final_face = [final_face; qiu_d];
	 
end

3.3 Mixed Noise Level

Function: Combine the planar fitting noise level and the minimum bounding sphere noise level with weights to obtain the mixed noise level.
Implementation:
matlab

hunhe_d = 0.4 * ping_d + 0.6 * final_face;

4.Noise Point Identification

Function: Identify the noise points under planar fitting, minimum bounding sphere fitting, and mixed fitting respectively according to the set noise probability thresholds.
Implementation:
matlab

pingyuzhi = 0.68;

qiuyuzhi = 0.38;

Hyuzhi = 0.50;

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

5.Visualization

Function: Visualize the simulated noise points and the original point cloud, as well as the noise points identified by different methods and the original point cloud in different windows respectively.
Implementation:
matlab

pcshowpair(ptCloudA, ptCloud); 

figure;

ptCloudss = pointCloud(ping_p);

pcshowpair(ptCloudss, ptCloud); 

figure;

ptClouds = pointCloud(qiu_p); 

pcshowpair(ptClouds, ptCloud); 

figure;

ptCloudsh = pointCloud(h_p);

pcshowpair(ptCloudsh, ptCloud);

Part3 Fig8_ProMeanShift_E

1.Introduction

This MATLAB script is designed to process point cloud data, add simulated Gaussian noise, and then use the mean shift algorithm to smooth the noisy point cloud. The following sections provide detailed information about the script, including dependencies, key algorithms, and implementation details.

2.Dependencies and Requirements

2.1 Software Requirements

MATLAB R2020b or later. The script uses functions such as uigetfile, importdata, pointCloud, knnsearch, normrnd, meanshift, etc., which are available in the standard MATLAB environment or require additional toolboxes.
Some custom functions (funzhong_zaoshengdu, Funpca, findNearestNeighbors, select) are used in the script. These functions should be available in the MATLAB path.

2.2 Input Data

The script supports point cloud data files in .asc, .txt, and .xlsx formats. The user can select a point cloud data file through a file dialog.

3.Key Algorithms and Implementation

3.1 Point Cloud Data Loading

matlab

[FileName, PathName] = uigetfile({'*.asc'; '*.txt'; '*.xlsx'}, 'Select point cloud data');

str = [PathName FileName];

pointdatas = importdata(str);

ptCloud = pointCloud(pointdatas);

The uigetfile function is used to open a file dialog, allowing the user to select a point cloud data file.

The importdata function reads the selected file and stores the data in pointdatas.

A pointCloud object is created using the imported data.

3.2 Simulated Noise Addition

matlab

a = ptCloud.Location;

sigma = 0.001;

s_noise = normrnd(0, sigma, round(0.2 * size(a, 1)), 3);

s_simulated = a(1:5:5 * round(0.2 * size(a, 1)), 1:3) + s_noise;

ss = [s_simulated; a(:, 1:3)];

ptCloudA = pointCloud(s_simulated);

ptCloudB = pointCloud(ss);

Gaussian noise with a mean of 0 and a standard deviation of sigma is generated.

A part of the original point cloud is selected every 5 points, and the generated noise is added to these points.
Two new pointCloud objects are created: ptCloudA contains only the noisy points, and ptCloudB contains both the original point cloud and the noisy points.

3.3 K - Nearest Neighbor Search

matlab

P = ss';

k = 25;

neighbors = transpose(knnsearch(transpose(P), transpose(P), 'k', k + 1));

The knnsearch function is used to find the k nearest neighbors of each point in the point cloud. The result is transposed and stored in neighbors.

3.4 Comprehensive Noise Degree Calculation

matlab

[hunhe_d, qiuzaosheng, bianfen, dist_m, dist, pn, pw] = funzhong_zaoshengdu(P, k, neighbors);

nf = pn';

The custom function funzhong_zaoshengdu is used to calculate the comprehensive noise degree and related parameters.

3.5 Principal Component Analysis (PCA) for Normal Vector Estimation

matlab

[pn, ~] = Funpca(P, k, neighbors);

nf = pn';

The custom function Funpca is used to perform principal component analysis on the data to estimate the normal vector of each point.

3.6 Mean Shift Clustering and Smoothing

matlab

sigmac = 0.3;

tol = 0.004;

[numPts, ~] = size(dataPts);

z = zeros(numPts, 3);

for n = 1:numPts
 
		% ...
 
		while error > tol
   
				% ...
 
				newx = newp * plo.Location;
 
				x_m = newx;

				x_f = x + hunhe_d(n).*dot(np,(x_m - x))*np;
 
				error = norm(newx - oldx);
		
    end
  
		zz(n, :) = x_f;
	
end
The mean shift algorithm is used to cluster the normal vectors of each point and its neighbors.
Gaussian kernel functions are used to calculate the weights based on direction and position differences.
The position of each point is iteratively updated until the error is less than the set threshold.

3.7 Point Cloud Visualization

matlab

pcshow(ptCloud);

pcshow(ptCloudB);

ptt2 = pointCloud(zz);

pcshow(ptt2);

The pcshow function is used to display the original point cloud, the point cloud with added noise, and the point cloud after mean shift smoothing.









