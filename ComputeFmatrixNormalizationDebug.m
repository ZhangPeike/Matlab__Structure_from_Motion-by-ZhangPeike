function [ F,Diag,SingularVec1,SingularVec2 ] = ComputeFmatrixNormalizationDebug(matches)
%��������ͼ������������
%���ߣ������
%���裺1.����SIFT����������ȡ��ƥ��
%2.ƥ���Լ���������Af=0
%3.��A��������ֵ�ֽ⣬SVD����С����ֵ��Ӧ������������f
%4.fת��ΪF��3*3����
%Version:2
%������޸ģ������˹�һ��
%matchesΪƥ���ԣ�4*n��n�Ե㣬һ��Ϊһ��
%F�Ĺ�һ����t�ĳ߶ȵ�Ӱ��
[~,b]=size(matches);
%���췽��
%2016��1��26��11:06:37����MyNormalizing����
%%
x1=matches(1:2,:);
x2=matches(3:4,:);
%��һ������һ�����������⣬���������ת����ĵ�
[T1,x1]=MyNormalizing(x1);
[T2,x2]=MyNormalizing(x2);
%%
%�����µ���ʽ����
%2016��2��7��11:49:31
%
%{
 A = [x2(1,:)'.*x1(1,:)'   x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
         x2(2,:)'.*x1(1,:)'   x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
         x1(1,:)'             x1(2,:)'            ones(npts,1) ];       

%}

A=[x2(1,1)*x1(1,1),x2(1,1)*x1(2,1),x2(1,1),x2(2,1)*x1(1,1),x2(2,1)*x1(2,1),x2(2,1),x1(1,1),x1(2,1),1];
for i=1:b-1
    A=[A;x2(1,1+i)*x1(1,1+i),x2(1,1+i)*x1(2,1+i),x2(1,1+i),x2(2,1+i)*x1(1,1+i),x2(2,1+i)*x1(2,1+i),x2(2,1+i),x1(1,1+i),x1(2,1+i),1];
end
%��������ֵ�ֽ�
Diag=svd(A);

[~,~,V]=svd(A);
SingularVec1=V(:,9);
SingularVec2=V(:,8);
%V�����һ�м�f��˳��9��Ԫ��
tempF = reshape(V(:,9),3,3)';
%{
tempF=[V(1,9),V(2,9),V(3,9);
   V(4,9),V(5,9),V(6,9);
   V(7,9),V(8,9),V(9,9)];
%}
%F��Լ��������rank(F)=2
[U1,S1,V1]=svd(tempF);
tempF=U1*diag([S1(1,1),S1(2,2),0])*V1';
F=T2'*tempF*T1;