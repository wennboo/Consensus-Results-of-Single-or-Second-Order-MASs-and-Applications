
%%Multiplicative Noisy Measurements
%%图的拓扑说明：所有的有向图（不一定是平衡图）的并集含生成树，此仿真含三个拓扑。
%%%周期为3，每一个周期内拓扑切换A0-A1-A2，联合含生成树
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
M=1;%每个拓扑下的迭代次数，固定为M
N=8000;%切换次数
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
%h=2;%周期为2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%拓扑图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%有公共的左特征向量
%%%有向连接1-2-3-4,1-3-4-1
A1=[0 0 0 1;
    2 0 0 0;
    1 1 0 0;
    0 1 2 0];
A1=0.2*A1;
%%%有向连接2-3-4-1,1-4,4-2
A2=[0 0 0 1;
    0 0 0 1;
    0 1 0 0;
    3 0 1 0];
A2=0.2*A2;
%%%有向连接1-2-3-4，3-2-4,4-1
A3=[0 0 0 1;
    3 0 1 0;
    0 2 0 0;
    0 2 1 0];
A3=0.2*A3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(:,:,1)=A1;
A(:,:,2)=A2;
A(:,:,3)=A3;
n=size(A);
n1=n(1);
n2=n(3);
E=-eye(n1);
I=eye(n1);
E1=eye(n1-1);
d=ones(n1,1);
E(:,1)=d%%%转换矩阵E1
Gama=E(2:n,:);
pai=[3 1 1 1]';
pai=pai/6;
v_end=pai'*v0
E=[pai';Gama]
inv(E)
S=inv(E)'*diag(pai)*inv(E)
S=S(2:n,2:n)
eig(S)
%%%%%%%%
%%beta求解
beta=[];
g2=[];
g3=[];
for i=1:1:n2
    Tem=A(:,:,i).^2;
    beta=[beta sum(Tem(:))];
    L(:,:,i)=(diag(sum(A(:,:,i),2))-A(:,:,i));
    L0(:,:,i)=inv(E)'*diag(pai)*L(:,:,i)*inv(E);
    H=L0(:,:,i);
    F(:,:,i)=H(2:n,2:n);
    H=eig(F(:,:,i)'+F(:,:,i));
    g2=[g2 min(H)];
    H=eig(F(:,:,i)'*inv(S)*F(:,:,i));
    g3=[g3 max(H)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%mode1的增益求解
%%%%%%%%%%%%%%%
alpha=diag(Gama'*S*Gama);
alpha=max(alpha)
beta1=max(beta);
delta1=alpha*beta1*sigma01^2*sigma11^2
delta2=alpha*beta1*sigma02^2*sigma12^2
g1=max(eig(S));
g2=min(g2);
g3=max(g3);
c21=(-k*g1*(delta1+g3)+sqrt(k^2*g1^2*(delta1+g3)^2+8*k*g1*g3*g2^2))/(4*g2*g3)
c22=(k-3)*g2/(k*(delta2+2*g3))
cmax=min(c21,c22)
c1max=g2*cmax^2/(k*g1)
c2=cmax*0.95
c1=g2*c2^2/(k*g1)

% L1=L(:,:,1)
% L1'*pai
% H1=[I I;-c1*L1 I-c2*L1];
% H1^N*yita0
%H1=[I I;-c1*L10 I-c2*L10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%figure
flag=1;
if(flag)
c=[c1 c2];
K=size(yita0);
yita00=yita0;
K=K(2)-1;
yita=[];
for i=1:N
    j=unidrnd(3);
    switch(j)
        case 1,
            yita1=fix_sta_sol(yita00,L(:,:,1),A(:,:,1),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
        case 2,
            yita1=fix_sta_sol(yita00,L(:,:,2),A(:,:,2),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
        case 3,
            yita1=fix_sta_sol(yita00,L(:,:,3),A(:,:,3),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
    end
end

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
x_end=x(:,end)
v_end=v(:,end)
end