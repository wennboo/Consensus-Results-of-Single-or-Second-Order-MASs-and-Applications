
%����ϵͳ״ֵ̬
%AΪ�ڽӾ���LΪ���Ͼ���deltaΪ������׼�sigmaΪ|fji(xj-xi)|<=sigma|xj-xi|
%%%%%%%�˳���fji(xj-xi)=abs(sigma��xj-xi��)
function [x,Su]=dis_rsd_sta_sol0(x0,L,A,M,c,delta,ave,sigma,K)
x=[];
n=length(x0);
E=eye(n);
D=[];
Su=[];
for j=1:n;
D=blkdiag(D,A(j,:));
end
for i=K+1:K+M
    y=[];
    for m1=1:1:n
        for m2=1:1:n
            y=[y x0(m2)-x0(m1)];
        end
    end
    y=sigma*y;
    y=diag(y);
    y=abs(y);
    kesi=ave+sqrt(delta)*randn(n^2,1);
    x1=(E-c*L)*x0+c*D*y*kesi;
    Su=[Su E-c*L];
    x=[x x1];
    x0=x1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%