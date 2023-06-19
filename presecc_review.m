clc
clear
close all;

% 信号的点数
t = 0:0.1:50-0.1;
% 散射点个数
a = zeros(100,1);
% 进动频率
f0 = 20e9;
% 距离
Rt = 100000+100*t;
% 坐标
x = 0:0.2:2-0.2;
xm =  repmat(x, 10,1);
xm = (xm.');
xm = xm(:);
y = 0:0.2:2-0.2;
ym =  repmat(y, 10,1);
ym = ym(:);
% 字典矩阵
Q = zeros(length(t),length(a));
phit = zeros(length(t),1);
c = 3e8;
a(1) = 1.2+2i;a(50) = 10+5i;a(9) = 2;a(8)=3;
theta = 5*pi/180;
beta = 401*pi/180+0.2*pi/180*t;
% for i = 1:length(t)
    phit = acos(sin(theta).*sin(beta).*cos(2*pi*f0*(t)+75*pi/180)./cos(theta)^2 + cos(theta).*cos(beta) - sin(theta)^3.*sin(beta)/cos(theta)^2);
% end
s = zeros(length(t),1);
for i = 1:length(t)
    for j = 1:length(a)
        phi = 4*pi*f0*(Rt(i)+xm(j)*cos(phit(i))+ym(j)*sin(phit(i))/c);
        Q(i,j) = exp(-1*1j*phi);
    end
end
s = Q*a;
[x_hat,residual_start,residual_end] = omp(Q, s, 4);
% [x_hat] = omp(s, Q, 4);
figure;
plot(abs(x_hat));
k = find(x_hat>1)