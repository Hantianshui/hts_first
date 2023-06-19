% 画高斯函数图像
clc
clear;
close all;
x = -5:0.001:5;
y1 (:,1)= Gaussian(x,0,1)';
x1 = -5:0.001:5;
y1(:,2) = Gaussian(x1,0,2)';
% plot(x,y0,'b');
area(x1,y1)
legend('方差为1','方差为2')
title('高斯分布概率密度积分图');
xlabel('角度偏离值')
ylabel('概率密度')
