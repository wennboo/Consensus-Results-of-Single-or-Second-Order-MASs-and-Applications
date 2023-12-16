
%返回系统状态值
%A为邻接矩阵，L为拉氏矩阵，delta为噪声标准差，sigma为|fji(xj-xi)|<=sigma|xj-xi|
%%%%%%%此程序fji(xj-xi)=abs(sigma（xj-xi）)
function yita=fix_sta_sol(yita0,L,A,M,c,delta,ave,sigma,K)
yita=[];
n=length(yita0);
n1=n/2;
x0=yita0(1:n1);
v0=yita0(n1+1:n);
E=eye(n1);
Zero1=zeros(n1,1);
D=blkdiag(A(1,:),A(2,:),A(3,:),A(4,:));
for i=K+1:K+M
    y=[];
    z=[];
    for m1=1:1:n1
        for m2=1:1:n1
            y=[y x0(m2)-x0(m1)];
            z=[z v0(m2)-v0(m1)];
        end
    end
    y=sigma(1)*y;
    y=diag(y);
    y=abs(y);
    z=sigma(2)*z;
    z=diag(z);
    z=abs(z);
    kesi1=ave+sqrt(delta(1))*randn(n1^2,1);
    kesi2=ave+sqrt(delta(2))*randn(n1^2,1);
    H1=[E E;-c(1)*L E-c(2)*L];
    W1=[Zero1;c(1)*D*y*kesi1];
    W2=[Zero1;c(2)*D*z*kesi2];
    yita1=H1*yita0+W1+W2;
    yita=[yita yita1];
    yita0=yita1;
    x0=yita0(1:n1);
    v0=yita0(n1+1:n);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%