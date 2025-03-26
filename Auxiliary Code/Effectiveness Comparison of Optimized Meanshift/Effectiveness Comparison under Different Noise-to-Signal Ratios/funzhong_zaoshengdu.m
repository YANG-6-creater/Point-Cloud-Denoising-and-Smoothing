function [hunhe_d,qiuzaosheng,bianfen,dist_m,dist,n,w] = funzhong_zaoshengdu(p, k, neighbors)%p为输入的点云
if nargin < 2
    error('no bandwidth specified')
end
if nargin < 3
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));%neighbor为一个索引矩阵，第一行代表第几个点，后8行代表K近邻的点。记录每个点及其周围的8个点
end
m = size(p,2);
%返回第2维的维度
n = zeros(3,m); %存放法线的矩阵
w = zeros(1,m);
dist=zeros(1,m);
dist_m=zeros(1,m);
bianfen = zeros(1,m);
qiuzaosheng=zeros(1,m);%存储基于最小球的噪声度
for i = 1:m
    x = p(:,neighbors(2:end, i));%x为个邻域点    ，3x8的矩阵
    p_bar=mean(x,2);%每一行求均值(一共三行)
    
    point=p(:,i);
%     r_qiu=sqrt(transpose(point-p_bar)*(point-p_bar));
%     r_qius=sqrt(transpose(x-repmat(p_bar,1,k))*(x-repmat(p_bar,1,k)));
    r_qiu=pdist2(point',p_bar');
    r_qius=pdist2(x',p_bar');
    qiuzaosheng(:,i)=r_qiu/(2*mean(r_qius)/sqrt(k)+r_qiu);
%     P =  (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k));%中心化样本矩阵，再计算协方差矩阵
%     P = 1/(k) * (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %邻域协方差矩阵P
    P=(x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k))./(size(x,2)-1);
    
    
    [V,D] = eig(P);%求P的特征值、特征向量。 D是对应的特征值对角矩阵，V是特征向量(因为协方差矩阵为实对称矩阵，故特征向量为单位正交向量)
    r=sqrt(diag(V));%
    bianfen(:,i)=r(1)/(r(1)+r(2)+r(3));
    
    [d0, idx] = min(diag(D)); %d0为最小特征值  idx为特征值的列数索引。diag():创建对角矩阵或获取矩阵的对角元素
 
    
    n(:,i) = V(:,idx);   % 最小特征值对应的特征向量为法矢，即法向量 
    dm=abs((x - repmat(p_bar,1,k))'*n(:,i))./sqrt(n(:,i)'*n(:,i));
    dist_m(:,i)=mean(abs((x - repmat(p_bar,1,k))'*n(:,i))./sqrt(n(:,i)'*n(:,i)));
    dist(:,i)=abs((p(:,i)-p_bar)'*n(:,i))/sqrt(n(:,i)'*n(:,i));
    
%     dist(:,i)=dist2plane(p(:,m),);
    %规定法矢方向指向
    flag = p(:,i) - p_bar;%由近邻点的平均点指向对应点的向量
    if dot(n(:,i),flag)<0%如果这个向量与法向量的数量积为负数（反向）
        n(:,i)=-n(:,i);%法向量取反向
    end
    if nargout > 1 
        w(1,i)=abs(d0)./sum(abs(diag(D)));%最小特征值的绝对值在协方差矩阵特征值绝对值的总和中占的比重
    end
end
pointdata=p';
[pointdata_rows,pointdata_list] = size(pointdata);
ping_d=dist./(dist+dist_m);%平面拟合噪声度 d/(d_m+d)
ping_d=ping_d';
% mean_dists = mean(dists); %求均值
% ss = std(dists,0,2); %求标准差
% max_dists = mean_dists+(1*ss); %定义距离
final_face = [];%最小球噪声度
ptCloudB = pointCloud(pointdata(:,1:3));
for ii = 1:pointdata_rows
    point = pointdata(ii,:);
    [indices,dists] = findNearestNeighbors(ptCloudB,point,k);
    p = select(ptCloudB,indices);
    p=p.Location;
    center=mean(p,1);
    center_dist=pdist2(point,center);
    center_dists=pdist2(point,p);
    qiu_d=center_dist/(2*mean(center_dists)/sqrt(k)+center_dist);
    final_face = [final_face;qiu_d];
%   indices_list = []; 
%     if dists(ii) < max_dists
%         final_face = [final_face;point];
%     end
end
hunhe_d=0.4*ping_d+0.6*final_face;