function [x,tim,Su]=con_rsd_sta_sol0(x0,t0,L,A,tao,c,delta,ave,sigma,K,deltat)
x=[];
tim=[];
t1=[];
Su=[];
n=length(x0);
E=eye(n);
D=[];
M=tao/deltat;
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
    x1=(E-deltat*c*L)*x0+sqrt(deltat)*c*D*y*kesi;
    Su=[Su E-deltat*c*L];
    t1=t0+deltat;
    x=[x x1];
    tim=[tim t1];
    x0=x1;
    t0=t1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%