% ����˹����ͼ��
clc
clear;
close all;
x = -5:0.001:5;
y1 (:,1)= Gaussian(x,0,1)';
x1 = -5:0.001:5;
y1(:,2) = Gaussian(x1,0,2)';
% plot(x,y0,'b');
area(x1,y1)
legend('����Ϊ1','����Ϊ2')
title('��˹�ֲ������ܶȻ���ͼ');
xlabel('�Ƕ�ƫ��ֵ')
ylabel('�����ܶ�')
