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
%%%%%每一个含生成树
%%%有向连接1-2-3-4,2-1
A1=[0 1 0 0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0];
A1=0.2*A1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%有公共的左特征向量
%%%有向连接1-3,1-2-4，4-1
A2=[0 0 0 1;
    1 0 0 0;
    1 0 0 0;
    0 1 0 0];
A2=0.2*A2;
%%%有向连接4-2-1,2-3-4
A3=[0 1 0 0;
    0 0 0 1;
    0 1 0 0;
    0 0 1 0];
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
E(:,1)=d;%%%转换矩阵E1
Gama=E(2:n,:);
beta=[];
bar_alpha1=[];
bar_alpha2=[];
bar_g1=[];
bar_g2=[];
c2max=[];
for i=1:1:n2
    Tem=A(:,:,i).^2;
    beta1=sum(Tem(:));
    beta=[beta beta1];
    L(:,:,i)=(diag(sum(A(:,:,i),2))-A(:,:,i));
    L0(:,:,i)=inv(E)*L(:,:,i)*E;
    H=L0(:,:,i);
    H=H(2:n,2:n);
    L1(:,:,i)=H;
    P=lyap(H',-E1);
    bar_P(:,:,i)=P;
    alpha=diag(Gama'*P*Gama);
    alpha=max(alpha);

    alpha1=alpha*beta1*sigma01^2*sigma11^2;
    bar_alpha1=[bar_alpha1 alpha1];
    alpha2=alpha*beta1*sigma02^2*sigma12^2;
    bar_alpha2=[bar_alpha2 alpha2];
    g1=max(eig(P));
    bar_g1=[bar_g1 g1]; 
    g2=max(eig(H'*P*H));
    bar_g2=[bar_g2 g2];
    c21=(-k*g1*(alpha1+g2)+sqrt(k^2*g1^2*(alpha1+g2)^2+8*k*g1*g2))/(4*g2);
    c22=(k-3)/(k*(alpha2+2*g2));
    c2max1=min(c21,c22);
    c2max=[c2max c2max1];
end
% beta
% bar_alpha
c2max=min(c2max)
c1max=c2max^2/(k*max(bar_g1))
% bar_g1
% bar_g2
c2=c2max*0.95
c1=c2^2/(k*max(bar_g1))
gamap=[];
mu=[];
mu1=[];
for i=1:1:n2
    P=bar_P(:,:,i);
    P1=[c1*E1 c1*P/c2;c1*P/c2 P];
    a11=bar_alpha1(i);
    a12=bar_alpha2(i);
    L11=L1(:,:,i);
    Omega11=(-c1^2/c2+a11*c1^2)*E1+c1^3/c2*L11'*P*L11;
    Omega12=-c1^2/c2*L11'*P+c1*c2*L11'*P*L11;
    Omega22=(-c2+a12*c2^2)*E1+2*c1/c2*P+c2^2*L11'*P*L11;
    Omega=[Omega11 Omega12;Omega12' Omega22];
    eig(Omega)
    H=eig(inv(P1)*Omega);
    gama1=1+max(H);
    gamap=[gamap gama1];
    for j=1:1:n2
        if(j~=i)
            Q=bar_P(:,:,j);
            Q1=[c1*E1 c1*Q/c2;c1*Q/c2 Q];
            mu2=max(eig(inv(Q1)*P1));
            mu1=[mu1 mu2];
        end
    end
    mu2=max(mu1);
    mu1=[];
    mu=[mu mu2];

end
gamap
mu
tao=[];
T=[];
for i=1:1:n2
    if(mu(i)>1)
        tao1=-log(mu(i))/log(gamap(i));
        tao=[tao tao1];
        T1=ceil(tao1);
        T=[T T1];
    else
        tao1=1;
        tao=[tao tao1];
        T=[T tao1];
    end
end
tao
T
flag=1;
if(flag)
c=[c1 c2];
K=size(yita0);
yita00=yita0;
K=K(2)-1;
yita=[];
M2=2
for i=1:M2
    for j=1:n2
    for i0=1:T(j)
            yita1=fix_sta_sol(yita00,L(:,:,j),A(:,:,j),M,c,sigma1,ave,sigma0,K);
            yita=[yita yita1];
            K=size(yita);
            yita00=yita(:,K(2));
            K=K(2)-1;
    end
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
length(x)