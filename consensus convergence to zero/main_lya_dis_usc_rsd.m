
%%Multiplicative Noisy Measurements
%%图的拓扑说明：所有的有向图（不一定是平衡图）的并集含生成树，此仿真含2个拓扑。
%%%周期为2，每一个周期内拓扑切换A0-A1，联合含生成树
clc;
clear;
close all;
%%%%

%%系统赋初值
%x0=[4,3,-3,-2]';
x0=[1.5,1,-0.5,-1]';
n0=length(x0);
E=eye(n0);
%%%%迭代次数
M=1;%每个拓扑下的迭代次数，固定为M
N=4000;

pai_inf=[0.3 0.2 0.3 0.2]';
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%噪声的设置
sigma1=0;%相对噪声强度|f(x)-f(y)|=mu|x-y|
delta=0;%用白噪声方差为1，用点乘实现
ave=0;
h=2;%周期为3
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%给出c验证gama的值是否属于(0,1)
%c=0.005;
%c=0.005;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%拓扑图
%%%G0与G1组成union-mode 1
%%%有向连接2-1；4-1
%%%有向连接2-1；4-1
A1=[0 0 1 0;
    0 0 0 0;
    0 0 0 0;
    0 1 0 0];

%%%有向连接4-3,4-2,1-3;
A2=[0 0 0 0;
    1 0 0 0;
    0 0 0 1;
    1 0 0 0];

A3=[0 1 1 0;
    0 0 0 0;
    0 1 0 0;
    0 0 0 0];
%A1=0.4*A1;
Au(:,:,1)=A1;
Au(:,:,2)=A2;
Au(:,:,3)=A3
%%二分法确定一个值尽可能大的c
c=0.01
n=size(Au);
n1=n(3);
%%flag控制后面运行
for i=1:1:n1
    Lu(:,:,i)=(diag(sum(Au(:,:,i),2))-Au(:,:,i));
end

%%flag控制后面运行
flag=1;
if(flag)
x=x0;
x00=x0;
K=size(x);
K=K(2)-1;
n0=length(x0);
swi_j=[];
Su=[];
for i=1:N
    j=unidrnd(2);
    swi_j=[swi_j j*ones(1,2*M)];
    switch(j)
        case 1,
            [x1,Su1]=dis_rsd_sta_sol0(x00,Lu(:,:,1),Au(:,:,1),M,c,delta,ave,sigma1,K);
            Su=[Su Su1];
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1;
            [x1,Su1]=dis_rsd_sta_sol0(x00,Lu(:,:,2),Au(:,:,2),M,c,delta,ave,sigma1,K);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1;
        case 2,
            [x1,Su1]=dis_rsd_sta_sol0(x00,Lu(:,:,2),Au(:,:,2),M,c,delta,ave,sigma1,K);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1;
            [x1,Su1]=dis_rsd_sta_sol0(x00,Lu(:,:,3),Au(:,:,3),M,c,delta,ave,sigma1,K);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1;
    end
end
n=length(x)
%title('Average path');  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%最后一次样本值仿真    
N1=length(x0);
N=length(Su);
N=floor(N/N1);

 Pai=[];
 pai1=pai_inf';
for i=1:N
    tem1=4*(N-i+1);
    tem0=4*(N-i)+1;
    pai0 =pai1*Su(:,tem0:tem1);
    Pai=[pai0' Pai];
    pai1=pai0; 
end


n=length(Pai);

vt=[];

for i=1:n 
    Qk=diag(Pai(:,i))-Pai(:,i)*Pai(:,i)';
    vt0=x(:,i)'*Qk*x(:,i);
    vt=[vt vt0];
end
save dis_lya vt x;



for i=1:N1
    plot(x(i,1:n-1),'k');
    hold on;
end
figure;
plot(vt(1:n),'r');
figure
plot(log(vt(1:n)),'r');
end



