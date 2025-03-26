function [clustCent,data2cluster,numClust] = meanshift(dataPts,bandWidth)
% �����޺�Ư�ƾ����㷨
%perform MeanShift Clustering of data using a flat kernel
%stopThresh = 1e-3*bandWidth; %������ֹ�ķ�ֵ
% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% bandWidth         - is bandwidth parameter (scalar)���Ǵ��������������
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
    error('�����������');
end


[numDim,numPts] = size(dataPts);%[ά�ȣ�����][����������]����,��ʼ��δ�������
% ��ʼ��
numInitPts = numPts; %��ʼ����δ�������
Initidx = 1:numPts; %��ʼ����δ���������
% clusterVotes = zeros(100,numPts,'uint16'); %��ʼ���Ĵ�ͶƱ
clusterVotes = zeros(5,numPts,'uint16'); %��ʼ���Ĵ�ͶƱ
beenVisitedFlag = zeros(1,numPts,'uint8'); %�ж�ѭ����ֹ����
numClust = 0; %��ʼ���Ĵ�����
% clustCent = zeros(numDim,100); %��ʼ���Ĵ�����
clustCent = zeros(numDim,5); %��ʼ���Ĵ�����
while numInitPts
    tempInd = ceil( (numInitPts-1e-6)*rand); %ѡ��������ӵ㣬��Ϊ��ʼ�������� %�˴����þ��ȳ�ֵ�Ż�
    stInd = Initidx(tempInd); %�õ���ԭʼ�����е�����
    myMean = dataPts(:,stInd); %�õ�����
    stopThresh = 1e-3*bandWidth; %������ֹ�ķ�ֵ
    thisClusterVotes = zeros(1,numPts,'uint16'); %����Ըôص�ͶƱͳ��
    while 1
        %�����ֵƯ������
        inInds = sum((repmat(myMean,1,numPts) - dataPts).^2) < bandWidth.^2; %�ڵ���߼����� %�˴�����k-d tree����
        thisClusterVotes(inInds) = thisClusterVotes(inInds)+1; %���ڵ�Ϊ�ô�����һ��ͶƱ
        myOldMean = myMean;
        myMean  = mean(dataPts(:,inInds),2); %�����µľ�����������
        beenVisitedFlag(inInds) = 1; %��ǵ����и���Ĳ�ѯ���
        if norm(myMean-myOldMean) < stopThresh %ѭ����ֹ����
            mergeWith = 0; %��ʼ���ϲ��صĲ���
            for cN = 1:numClust
                distToOther = norm(myMean-clustCent(:,cN)); %�����������ĵľ���
                if distToOther < bandWidth/2
                    mergeWith = cN;
                    break
                end
            end
            if mergeWith > 0
                clustCent(:,mergeWith) = 0.5*(myMean+clustCent(:,mergeWith)); %���յ�clustCent�����¼�������
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