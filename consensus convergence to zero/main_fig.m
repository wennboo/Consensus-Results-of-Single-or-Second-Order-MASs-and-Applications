
%%Multiplicative Noisy Measurements
%%图的拓扑说明：所有的有向图（不一定是平衡图）的并集含生成树，此仿真含2个拓扑。
%%%周期为2，每一个周期内拓扑切换A0-A1，联合含生成树
clc;
clear;
close all;
%%%%

val_struct=load('dis_lya');
val_names=fieldnames(val_struct)
dis_x=getfield(val_struct,val_names{1});
dis_vt=getfield(val_struct,val_names{2});

val_struct=load('con_lya');
val_names=fieldnames(val_struct)
con_x=getfield(val_struct,val_names{1});
con_tim=getfield(val_struct,val_names{2});
con_vt=getfield(val_struct,val_names{3});


subplot(2,2,1)
plot(dis_vt,'r');
subplot(2,2,2)
plot(log(dis_vt),'b');


[m0,n0]=size(con_x)
[m1,n1]=size(con_vt)
n=min(n0,n1);
subplot(2,2,3)
plot(con_tim(1:n),con_vt(1:n),'r');
subplot(2,2,4)
plot(con_tim(1:n),log(con_vt(1:n)),'b');




