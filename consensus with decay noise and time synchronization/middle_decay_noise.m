
%%����ʽ˥�����Բ������������� �ȶ�������
%%ͼ������˵�������е�����ͼ�Ĳ��������������˷��溬2�����ˡ�
%%%����Ϊ2��ÿһ�������������л�A0-A1�����Ϻ�������
clc;
clear;
close all;
%%%%

%%ϵͳ����ֵ
%x0=[4,3,-3,-2]';
x0=[3,-1,-2,-4]';
N0=length(x0);
E=eye(N0);
%%%%��������
M=1;%ÿ�������µĵ����������̶�ΪM
K=1;
N=5000;
flag=1;
flag1=1;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ͼ
%%%G0��G1���union-mode 1
%%%��������2-1��4-1
%%%��������2-1��4-1
A1=[0 0 1 0;
    0 0 0 0;
    0 0 0 0;
    0 1 0 0];

%%%��������4-3,4-2,1-3;
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
%%flag���ƺ�������
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

