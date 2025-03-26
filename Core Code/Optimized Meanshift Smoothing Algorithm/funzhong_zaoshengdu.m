function [hunhe_d,qiuzaosheng,bianfen,dist_m,dist,n,w] = funzhong_zaoshengdu(p, k, neighbors)%pΪ����ĵ���
if nargin < 2
    error('no bandwidth specified')
end
if nargin < 3
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));%neighborΪһ���������󣬵�һ�д���ڼ����㣬��8�д���K���ڵĵ㡣��¼ÿ���㼰����Χ��8����
end
m = size(p,2);
%���ص�2ά��ά��
n = zeros(3,m); %��ŷ��ߵľ���
w = zeros(1,m);
dist=zeros(1,m);
dist_m=zeros(1,m);
bianfen = zeros(1,m);
qiuzaosheng=zeros(1,m);%�洢������С���������
for i = 1:m
    x = p(:,neighbors(2:end, i));%xΪ�������    ��3x8�ľ���
    p_bar=mean(x,2);%ÿһ�����ֵ(һ������)
    
    point=p(:,i);
%     r_qiu=sqrt(transpose(point-p_bar)*(point-p_bar));
%     r_qius=sqrt(transpose(x-repmat(p_bar,1,k))*(x-repmat(p_bar,1,k)));
    r_qiu=pdist2(point',p_bar');
    r_qius=pdist2(x',p_bar');
    qiuzaosheng(:,i)=r_qiu/(2*mean(r_qius)/sqrt(k)+r_qiu);
%     P =  (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k));%���Ļ����������ټ���Э�������
%     P = 1/(k) * (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %����Э�������P
    P=(x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k))./(size(x,2)-1);
    
    
    [V,D] = eig(P);%��P������ֵ������������ D�Ƕ�Ӧ������ֵ�ԽǾ���V����������(��ΪЭ�������Ϊʵ�Գƾ��󣬹���������Ϊ��λ��������)
    r=sqrt(diag(V));%
    bianfen(:,i)=r(1)/(r(1)+r(2)+r(3));
    
    [d0, idx] = min(diag(D)); %d0Ϊ��С����ֵ  idxΪ����ֵ������������diag():�����ԽǾ�����ȡ����ĶԽ�Ԫ��
 
    
    n(:,i) = V(:,idx);   % ��С����ֵ��Ӧ����������Ϊ��ʸ���������� 
    dm=abs((x - repmat(p_bar,1,k))'*n(:,i))./sqrt(n(:,i)'*n(:,i));
    dist_m(:,i)=mean(abs((x - repmat(p_bar,1,k))'*n(:,i))./sqrt(n(:,i)'*n(:,i)));
    dist(:,i)=abs((p(:,i)-p_bar)'*n(:,i))/sqrt(n(:,i)'*n(:,i));
    
%     dist(:,i)=dist2plane(p(:,m),);
    %�涨��ʸ����ָ��
    flag = p(:,i) - p_bar;%�ɽ��ڵ��ƽ����ָ���Ӧ�������
    if dot(n(:,i),flag)<0%�����������뷨������������Ϊ����������
        n(:,i)=-n(:,i);%������ȡ����
    end
    if nargout > 1 
        w(1,i)=abs(d0)./sum(abs(diag(D)));%��С����ֵ�ľ���ֵ��Э�����������ֵ����ֵ���ܺ���ռ�ı���
    end
end
pointdata=p';
[pointdata_rows,pointdata_list] = size(pointdata);
ping_d=dist./(dist+dist_m);%ƽ����������� d/(d_m+d)
ping_d=ping_d';
% mean_dists = mean(dists); %���ֵ
% ss = std(dists,0,2); %���׼��
% max_dists = mean_dists+(1*ss); %�������
final_face = [];%��С��������
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