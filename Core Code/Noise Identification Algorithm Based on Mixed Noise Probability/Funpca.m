% PCA��Ԫ������������
% ���룺
% p:3*n����ֵ����
% k:k���ڲ���
% neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));
% neighborsһ���ȱʡ����֮ǰ����k������ȡ����Ҳ��ֱ�ӵ��ã��������Ч��
% ���
% n:��ʸ���ѹ涨������������ϳ���ƽ��ָ���ѯ��
% w:�����������ʵĲ����������Mark P,et al. Multi-scale Feature Extraction on Point-Sampled Surfaces[J]. Computer Graphics Forum, 2010, 22(3): 281-289.
 
function [n,w] = Funpca(p, k, neighbors)%pΪ����ĵ���
if nargin < 2
    error('no bandwidth specified')
end
if nargin < 3
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));%neighborΪһ���������󣬵�һ�д���ڼ����㣬��8�д���K���ڵĵ㡣��¼ÿ���㼰����Χ��8����
end
m = size(p,2);%���ص�2ά��ά��
n = zeros(3,m); %��ŷ��ߵľ���
w = zeros(1,m);
for i = 1:m
    x = p(:,neighbors(2:end, i));%xΪ8�������    ��3x8�ľ���
    p_bar=mean(x,2);%ÿһ�����ֵ(һ������)
    
%     P =  (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k));%���Ļ����������ټ���Э�������
%     P = 1/(k) * (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %����Э�������P
    P=(x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k))./(size(x,2)-1);
    
    [V,D] = eig(P);%��P������ֵ������������ D�Ƕ�Ӧ������ֵ�ԽǾ���V����������(��ΪЭ�������Ϊʵ�Գƾ��󣬹���������Ϊ��λ��������)
    
    [d0, idx] = min(diag(D)); %d0Ϊ��С����ֵ  idxΪ����ֵ������������diag():�����ԽǾ�����ȡ����ĶԽ�Ԫ��
 
    
    n(:,i) = V(:,idx);   % ��С����ֵ��Ӧ����������Ϊ��ʸ����������    
    
    %�涨��ʸ����ָ��
    flag = p(:,i) - p_bar;%�ɽ��ڵ��ƽ����ָ���Ӧ�������
    if dot(n(:,i),flag)<0%�����������뷨������������Ϊ����������
        n(:,i)=-n(:,i);%������ȡ����
    end
    if nargout > 1 
        w(1,i)=abs(d0)./sum(abs(diag(D)));%��С����ֵ�ľ���ֵ��Э�����������ֵ����ֵ���ܺ���ռ�ı���,����
    end
end