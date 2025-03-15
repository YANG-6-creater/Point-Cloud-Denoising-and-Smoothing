function [clustCent,data2cluster,numClust] = meanshift(dataPts,bandWidth)
% 经典无核漂移聚类算法
%perform MeanShift Clustering of data using a flat kernel
%stopThresh = 1e-3*bandWidth; %迭代中止的阀值
% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% bandWidth         - is bandwidth parameter (scalar)：是带宽参数（标量）
% plotFlag          - display output if 2 or 3 D    (logical)
% ---OUTPUT---
% clustCent         - is locations of cluster centers (numDim x numClust)
% data2cluster      - for every data point which cluster it belongs to (numPts)
% cluster2dataCell  - for every cluster which points are in it (numClust)
% 
% Bryan Feldman 02/24/06
% MeanShift first appears in
% K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a
% Density Function, with Applications in Pattern Recognition"


%*** Check input ****
if nargin < 2
    error('输入参数不足');
end


[numDim,numPts] = size(dataPts);%[维度，点数][行数，列数]列数,初始化未聚类点数
% 初始化
numInitPts = numPts; %初始化的未聚类点数
Initidx = 1:numPts; %初始化的未聚类点索引
% clusterVotes = zeros(100,numPts,'uint16'); %初始化的簇投票
clusterVotes = zeros(5,numPts,'uint16'); %初始化的簇投票
beenVisitedFlag = zeros(1,numPts,'uint8'); %判断循环中止条件
numClust = 0; %初始化的簇总数
% clustCent = zeros(numDim,100); %初始化的簇中心
clustCent = zeros(numDim,5); %初始化的簇中心
while numInitPts
    tempInd = ceil( (numInitPts-1e-6)*rand); %选择随机种子点，作为初始聚类中心 %此处可用均匀初值优化
    stInd = Initidx(tempInd); %该点在原始数据中的索引
    myMean = dataPts(:,stInd); %该点坐标
    stopThresh = 1e-3*bandWidth; %迭代中止的阀值
    thisClusterVotes = zeros(1,numPts,'uint16'); %各点对该簇的投票统计
    while 1
        %计算均值漂移向量
        inInds = sum((repmat(myMean,1,numPts) - dataPts).^2) < bandWidth.^2; %内点的逻辑索引 %此处可由k-d tree加速
        thisClusterVotes(inInds) = thisClusterVotes(inInds)+1; %各内点为该簇增加一个投票
        myOldMean = myMean;
        myMean  = mean(dataPts(:,inInds),2); %计算新的聚类中心坐标
        beenVisitedFlag(inInds) = 1; %标记点云中各点的查询情况
        if norm(myMean-myOldMean) < stopThresh %循环中止条件
            mergeWith = 0; %初始化合并簇的参数
            for cN = 1:numClust
                distToOther = norm(myMean-clustCent(:,cN)); %与其它簇中心的距离
                if distToOther < bandWidth/2
                    mergeWith = cN;
                    break
                end
            end
            if mergeWith > 0
                clustCent(:,mergeWith) = 0.5*(myMean+clustCent(:,mergeWith)); %最终的clustCent会重新计算所以
                clusterVotes(mergeWith,:) = clusterVotes(mergeWith,:) + thisClusterVotes;
            else
                numClust = numClust+1;
                clustCent(:,numClust) = myMean;
                clusterVotes(numClust,:) = thisClusterVotes;
            end
            break
        end
    end
    Initidx = find(beenVisitedFlag == 0);
    numInitPts = length(Initidx);
end
% if numClust<100
if numClust<5
    clusterVotes(numClust+1:end,:) = [];
end
[~,data2cluster] = max(clusterVotes,[],1);
end