
%%多项式衰减乘性测量噪声主程序 稳定性情形
%%图的拓扑说明：所有的有向图的并集含生成树，此仿真含2个拓扑。
%%%周期为2，每一个周期内拓扑切换A0-A1，联合含生成树
clc;
clear;
close all;
%%%%

%%系统赋初值
%x0=[4,3,-3,-2]';
x0=[3,-1,-2,-4]';
N0=length(x0);
E=eye(N0);
%%%%迭代次数
M=1;%每个拓扑下的迭代次数，固定为M
K=1;
N=5000;
flag=1;
flag1=1;
%%%%%%%%%%%%%%%%%%%%%%%%


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
%%%%%%%%%%%%%%%%%%%%%%
bar_gam=1;
hat_gam=0.7;
%%%%%%%%%%%%%%%%%%%%%%%
n=size(Au);
n1=n(3);
%%flag控制后面运行
for i=1:1:n1
    Lu(:,:,i)=(diag(sum(Au(:,:,i),2))-Au(:,:,i));
end

x00=x0;
x=[];
err=[];
for k=1:N
    j=unidrnd(2);
    switch(j)
        case 1,
            x1=decay_sta_sol(x00,Lu(:,:,1),Au(:,:,1),bar_gam,hat_gam,K,M);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1
            x1=decay_sta_sol(x00,Lu(:,:,2),Au(:,:,2),bar_gam,hat_gam,K,M);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1
        case 2,
            x1=decay_sta_sol(x00,Lu(:,:,2),Au(:,:,2),bar_gam,hat_gam,K,M);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1
            x1=decay_sta_sol(x00,Lu(:,:,3),Au(:,:,3),bar_gam,hat_gam,K,M);
            x=[x x1];
            K=size(x);
            x00=x(:,K(2));
            K=K(2)-1
    end
    
   
                  
end

[N0,N]=size(x);
x(:,end-10:end)

for i=1:N
    tem=max(x(:,i))-min(x(:,i));
    err=[err tem];
end


for i=1:N0
    plot(x(i,1:N),'k');
    hold on;
end

figure;
plot(log(err),'g');

save ('middle_gama','x','err');

