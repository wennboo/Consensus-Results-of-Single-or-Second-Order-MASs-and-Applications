
%%Multiplicative Noisy Measurements
%%图的拓扑说明：所有的有向图（不一定是平衡图）的并集含生成树，此仿真含三个拓扑。
%%%周期为3，每一个周期内拓扑切换A0-A1-A2，联合含生成树
clc;
clear;
close all;
flag=1;
if(flag)
    format long;
else
    format short;
end

%%系统赋初值
%x0=[4,3,-3,-2]';
x0=[2,1,-1,-2]';
v0=[0.5,0.2,0.1,-0.2]';
yita0=[x0;v0];
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
%%%有向连接1-2-3-4,3-1,2-1
A1=[0 0 1 0;
    1 0 0 0;
    1 0 0 0;
    1 0 0 0];
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
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%转换矩阵E~
E(:,1)=d;
Gama=E(2:n,:);
pai=[3 1 1 1]';
pai=pai/6;
%%%%%转换矩阵E-
bar_E=[pai';Gama];
%%%%%S求解
S=inv(bar_E)'*diag(pai)*inv(bar_E);
S=S(2:n,2:n);
%%%%%%%%
%%beta求解
beta=[];
g2=[];
g3=[];
for i=1:1:n2
    Tem=A(:,:,i).^2;
    beta=[beta sum(Tem(:))];
    L(:,:,i)=(diag(sum(A(:,:,i),2))-A(:,:,i));
    %%%%%%%%%%%%%%%%%%%%
    %%spanning tree
    if(i==1)
    L11=inv(E)*L(:,:,i)*E;
    L11=L11(2:n,2:n);
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %%%strongly connected with common pai
    if(i>1)
    L0(:,:,i-1)=inv(bar_E)'*diag(pai)*L(:,:,i)*inv(bar_E);
    H=L0(:,:,i-1);
    F(:,:,i-1)=H(2:n,2:n);
    H=eig(F(:,:,i-1)'+F(:,:,i-1));
    g2=[g2 min(H)];
    H=eig(F(:,:,i-1)'*inv(S)*F(:,:,i-1));
    g3=[g3 max(H)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%mode1的增益求解
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
P=lyap(L11',-E1);
a1=diag(Gama'*P*Gama);
a1=max(a1);
a11=a1*beta(1)*sigma01^2*sigma11^2;
a12=a1*beta(1)*sigma02^2*sigma12^2;
g11=max(eig(P));
g12=max(eig(L11'*P*L11));
c121=(-k*g11*(a11+g12)+sqrt(k^2*g11^2*(a11+g12)^2+8*k*g11*g12))/(4*g12);
c122=(k-3)/(k*(a12+2*g12));
cmax1=min(c121,c122);
%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=diag(Gama'*S*Gama);
alpha=max(alpha);
beta1=max(beta(2:n2));
delta1=alpha*beta1*sigma01^2*sigma11^2;
delta2=alpha*beta1*sigma02^2*sigma12^2;
g1=max(eig(S));
g2=min(g2);
g3=max(g3);
c211=(-k*g1*(delta1+g3)+sqrt(k^2*g1^2*(delta1+g3)^2+8*k*g1*g3*g2^2))/(4*g2*g3);
c222=(k-3)*g2/(k*(delta2+2*g3));
cmax2=min(c211,c222);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%仿真增益选择
cmax=min(cmax1,cmax2);
c2=cmax*0.95
c11=c2^2/(k*g11);
c21=g2*c2^2/(k*g1);
c1=min(c11,c21)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%正定矩阵P~求解
P1=[c1*E1 c1*P/c2;c1*P/c2 P]
eig(P1)
%P2=[c1*g2*E1 c1*S/c2;c1*S/c2 S]
P2=[c1*g2*E1 c1*S/c2;c1*S/c2 S]
eig(P2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%tao1与tao2的求解
gamap=[];
for i=1:1:n2
    if(i==1)
        Omega11=(-c1^2/c2+a11*c1^2)*E1+c1^3/c2*L11'*P*L11;
        Omega12=-c1^2/c2*L11'*P+c1*c2*L11'*P*L11;
        Omega22=(-c2+a12*c2^2)*E1+2*c1/c2*P+c2^2*L11'*P*L11;
        Omega=[Omega11 Omega12;Omega12' Omega22];
        eig(Omega)
        H=eig(inv(P1)*Omega);
        gama1=1+max(H);
        gamap=[gamap gama1];
    end
    if(i>1)
        alpha1=alpha*beta(i)*sigma01^2*sigma11^2;
        alpha2=alpha*beta(i)*sigma02^2*sigma12^2;
        H=F(:,:,i-1)'+F(:,:,i-1);
        H2=F(:,:,i-1)'*inv(S)*F(:,:,i-1);
        Omega11=-c1^2/c2*H+alpha1*c1^2*E1+c1^3/c2*H2;
        Omega12=c1*g2*E1-c1*H-c1^2/c2*F(:,:,i-1)'+c1*c2*H2;
        Omega22=2*c1/c2*S+(c1*g2+alpha2*c2^2)*E1-(c1+c2)*H+c2^2*H2;
        Omega=[Omega11 Omega12;Omega12' Omega22];
        eig(Omega11)
        eig(Omega22)
        eig(Omega)
        H=eig(inv(P2)*Omega);
        gama1=1+max(H);
        gamap=[gamap gama1];
    end
end
gamap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%A1为mode1，A2-A4构成mode2；
gama1=gamap(1)
gama2=max(gamap(2:n2))

% L21=inv(E)*L(:,:,2)*E;
% L21=L21(2:n,2:n);
% P=lyap(L21',-E1);
% P2=[c1*E1 c1*P/c2;c1*P/c2 P]
Q2=P2
eig(Q2)

Q1=P1
eig(Q1)

mu1=max(eig(inv(Q2)*Q1))
mu2=max(eig(inv(Q1)*Q2))
if(mu1>1)
    tao1=-log(mu1)/log(gama1);
else
    tao1=1;
end
if(mu2>1)
    tao2=-log(mu2)/log(gama2);
else
    tao2=1;
end
tao1
tao2
T1=ceil(tao1)
T2=ceil(tao2)
    
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
M2=1;
for i=1:M2
    for i0=1:T1
            yita1=fix_sta_sol(yita00,L(:,:,1),A(:,:,1),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
    end
    for i1=1:T2
           j=unidrnd(3);
        switch(j)
        case 1,
            yita1=fix_sta_sol(yita00,L(:,:,2),A(:,:,2),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
        case 2,
            yita1=fix_sta_sol(yita00,L(:,:,3),A(:,:,3),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
        case 3,
            yita1=fix_sta_sol(yita00,L(:,:,4),A(:,:,4),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
        end
    end
end

x=yita(1:n,:);
v=yita(n+1:2*n,:);
for i=1:n
    plot(x(i,:),'k');
    hold on;
end
figure
for i=1:n
    plot(v(i,:),'k');
    hold on;
end
x_end=x(:,end)
v_end=v(:,end)
end