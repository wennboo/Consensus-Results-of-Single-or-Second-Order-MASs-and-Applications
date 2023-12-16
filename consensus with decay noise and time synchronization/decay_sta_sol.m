
%返回系统状态值
%A为邻接矩阵，L为拉氏矩阵，delta为噪声标准差，sigma为|fji(xj-xi)|<=sigma|xj-xi|
%%%%%%%此程序fji(xj-xi)=abs(sigma（xj-xi）)
function x=decay_sta_sol(x0,L,A,bar_gam,hat_gam,K,M)
x=[];
n=length(x0);
E=eye(n);
for i=K+1:K+M
    ak=(1+i)^(-hat_gam);
    Delta=ak*(1+i)^(-bar_gam)*A;
    x1=(E-ak*L+Delta)*x0;
    x=[x x1];
    x0=x1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%