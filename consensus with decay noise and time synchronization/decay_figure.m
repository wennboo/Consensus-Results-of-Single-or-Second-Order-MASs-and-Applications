clc;
clear;
close all;
flag=1;
if(flag)
    format long;
else
    format short;
end

%%%%%%%%%%%%%%%%%
load ('middle_gama','x','err');
[N0,N]=size(x);
for i=1:N0
    plot(x(i,1:N),'k');
    hold on;
end
figure
load ('small_gama','x','err');
plot(log(err),'g');
hold on;
load ('middle_gama','x','err');
plot(log(err),'r');
hold on;
load ('big_gama','x','err');
plot(log(err),'b');
