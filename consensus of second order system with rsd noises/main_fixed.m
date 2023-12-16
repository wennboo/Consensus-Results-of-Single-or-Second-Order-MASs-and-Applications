%%Multiplicative Noisy Measurements for second order
%%ͼ������˵����fixed��
%%%
clc;
clear;
close all;
%%%%

%%ϵͳ����ֵ
%x0=[4,3,-3,-2]';
x0=[2,1,-1,-2]';
v0=[0.5,0.2,0.1,-0.2]';
yita0=[x0;v0]
n0=length(x0);
E=eye(n0);
%%%%��������
M=500;%ÿ�������µĵ����������̶�ΪM
% N=500;%�л�����
% M1=1;%1000����


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%����������
sigma01=0.5;%�������ǿ��|f(x)-f(y)|=mu|x-y|
sigma02=0.5;
sigma0=[sigma01 sigma02];
sigma11=2;%�ð���������Ϊ1���õ��ʵ��
sigma12=2;
sigma1=[sigma11 sigma12];
ave=0;%%%%������ֵ
k=4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ͼ
%%%mode 1
%%%��������1-2-3-4,3-1,2-1
A1=[0 1 0 0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0];
A1=0.2*A1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�й���������������
%%%��������1-2-3-4-1��1-3,2-4
A2=[0 0 0 1;
    2 0 0 0;
    1 1 0 0;
    0 1 2 0];
A2=0.2*A2;
%%%��������2-4-3-2,1-4,4-1
A3=[0 0 0 1;
    0 0 0 1;
    0 1 0 0;
    3 0 1 0];
A3=0.2*A3;
%%%��������1-4-3-1��
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
E(:,1)=d;%%%ת������E1
%%%%%%%%
%%beta���
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
%%%mode1���������
L10=L(:,:,1);%%%%%%%%mode1��L����
L11= L1(:,:,1);%%%%״̬�任���L~
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
