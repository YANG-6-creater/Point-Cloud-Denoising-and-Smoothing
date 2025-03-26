function [dist_m,dist,n,w] = fajuli(p, k, neighbors)%pΪ����ĵ���
if nargin < 2
    error('no bandwidth specified')
end
if nargin < 3
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));%neighborΪһ���������󣬵�һ�д���ڼ����㣬��8�д���K���ڵĵ㡣��¼ÿ���㼰����Χ��8����
end
m = size(p,2);%���ص�2ά��ά��
n = zeros(3,m); %��ŷ��ߵľ���
w = zeros(1,m);
dist=zeros(1,m);
dist_m=zeros(1,m);%
% bianfen = zeros(1,m);%������
for i = 1:m
    x = p(:,neighbors(2:end, i));%xΪ�������    ��3x8�ľ���
    p_bar=mean(x,2);%ÿһ�����ֵ(һ������)
    
%     P =  (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k));%���Ļ����������ټ���Э�������
%     P = 1/(k) * (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %����Э�������P
    P=(x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k))./(size(x,2)-1);
    
    [V,D] = eig(P);%��P������ֵ������������ D�Ƕ�Ӧ������ֵ�ԽǾ���V����������(��ΪЭ�������Ϊʵ�Գƾ��󣬹���������Ϊ��λ��������)
%     r=sqrt(diag(V));%q��������
%     bianfen(:,i)=r(1)/(r(1)+r(2)+r(3));
    
    [d0, idx] = min(diag(D)); %d0Ϊ��С����ֵ  idxΪ����ֵ������������diag():�����ԽǾ�����ȡ����ĶԽ�Ԫ��
    n(:,i) = V(:,idx);   % ��С����ֵ��Ӧ����������Ϊ��ʸ���������� 
    dm=abs((x - repmat(p_bar,1,k))'*n(:,i))./sqrt(n(:,i)'*n(:,i));
    dist_m(:,i)=mean(abs((x - repmat(p_bar,1,k))'*n(:,i))./sqrt(n(:,i)'*n(:,i)));
    dist(:,i)=abs((p(:,i)-p_bar)'*n(:,i))/sqrt(n(:,i)'*n(:,i));
    
%   dist(:,i)=dist2plane(p(:,m),);
    %�涨��ʸ����ָ��
    flag = p(:,i) - p_bar;%�ɽ��ڵ��ƽ����ָ���Ӧ�������
    if dot(n(:,i),flag)<0%�����������뷨������������Ϊ����������
        n(:,i)=-n(:,i);%������ȡ����
    end
    if nargout > 1 
        w(1,i)=abs(d0)./sum(abs(diag(D)));%��С����ֵ�ľ���ֵ��Э�����������ֵ����ֵ���ܺ���ռ�ı���
    end
end