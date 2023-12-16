%%Multiplicative Noisy Measurements for second order
%%图的拓扑说明：fixed。
%%%
clc;
clear;
close all;
%%%%

%%系统赋初值
%x0=[4,3,-3,-2]';
x0=[2,1,-1,-2]';
v0=[0.5,0.2,0.1,-0.2]';
yita0=[x0;v0]
n0=length(x0);
E=eye(n0);
%%%%迭代次数
M=500;%每个拓扑下的迭代次数，固定为M
% N=500;%切换次数
% M1=1;%1000采样


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%噪声的设置
sigma01=0.5;%相对噪声强度|f(x)-f(y)|=mu|x-y|
sigma02=0.5;
sigma0=[sigma01 sigma02];
sigma11=2;%用白噪声方差为1，用点乘实现
sigma12=2;
sigma1=[sigma11 sigma12];
ave=0;%%%%噪声均值
k=4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%拓扑图
%%%mode 1
%%%有向连接1-2-3-4,3-1,2-1
A1=[0 1 0 0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0];
A1=0.2*A1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%有公共的左特征向量
%%%有向连接1-2-3-4-1，1-3,2-4
A2=[0 0 0 1;
    2 0 0 0;
    1 1 0 0;
    0 1 2 0];
A2=0.2*A2;
%%%有向连接2-4-3-2,1-4,4-1
A3=[0 0 0 1;
    0 0 0 1;
    0 1 0 0;
    3 0 1 0];
A3=0.2*A3;
%%%有向连接1-4-3-1，
A4=[0 0 0 1;
    3 0 1 0;
    0 2 0 0;
    0 2 1 0];
A4=0.2*A4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(:,:,1)=A1;
A(:,:,2)=A2;
A(:,:,3)=A3;
A(:,:,4)=A4;
n=size(A);
n1=n(1);
n2=n(3);
E=-eye(n1);
I=eye(n1);
E1=eye(n1-1);
d=ones(n1,1);
E(:,1)=d;%%%转换矩阵E1
%%%%%%%%
%%beta求解
beta=[];
for i=1:1:n2
    Tem=A(:,:,i).^2;
    beta=[beta sum(Tem(:))];
    L(:,:,i)=(diag(sum(A(:,:,i),2))-A(:,:,i));
    L0(:,:,i)=inv(E)*L(:,:,i)*E;
    H=L0(:,:,i);
    L1(:,:,i)=H(2:n,2:n);
end
beta
Gama=E(2:n,:);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%mode1的增益求解
L10=L(:,:,1);%%%%%%%%mode1的L矩阵
L11= L1(:,:,1);%%%%状态变换后的L~
%%%%%%%%%%%%%%%
P=lyap(L11',-E1)
alpha=diag(Gama'*P*Gama);
alpha=max(alpha)
alpha1=alpha*beta(1)*sigma01^2*sigma11^2
alpha2=alpha*beta(1)*sigma02^2*sigma12^2
g1=max(eig(P))
g2=max(eig(L11'*P*L11))
c21=(-k*g1*(alpha1+g2)+sqrt(k^2*g1^2*(alpha1+g2)^2+8*k*g1*g2))/(4*g2)
c22=(k-3)/(k*(alpha2+2*g2))
cmax=min(c21,c22)
c1max=cmax^2/(k*g1)
c2=cmax*0.95
c1=c2^2/(k*g1)
%H1=[I I;-c1*L10 I-c2*L10];
c=[c1 c2];
K=size(yita0);
K=K(2)-1;
yita=fix_sta_sol(yita0,L10,A1,M,c,sigma1,ave,sigma0,K);
x=yita(1:n,:);
v=yita(n+1:2*n,:);
subplot(2,1,1);
for i=1:n
    plot(x(i,:),'k');
    hold on;
end
subplot(2,1,2);
for i=1:n
    plot(v(i,:),'k');
    hold on;
end
