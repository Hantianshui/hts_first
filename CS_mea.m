clc;
clear;
close;
%设置参数
c = 3e8;
f0 = 20e9;
fs = 2048;
Rt = 50e3;
%将二维区域划分为多个散射中心
rx = 0:0.06:3;
ry = 0:0.06:3;
%假设采样时间内的姿态角变化序列
deltaphi = 5/(0.25*fs);
phi = 50:deltaphi:55;
% phi = 50;
%构造散射中心模型回波（直接用坐标构建）
xm = [0,3,3,1];
ym = [1.5,3,0,1];
a = [5,3,3,1.5];
s1 = zeros(1,0.25*fs);
s = zeros(1,0.25*fs);
for i = 1:0.25*fs
    for k = 1:size(a,2)
        s1(i) = s1(i)+a(k).*(exp(-1j*4*pi*f0*((Rt+xm(k)*cosd(phi(i))+ym(k)*sind(phi(i)))/c)));
    end
end
Q = zeros(0.25*fs,size(rx,2)*size(ry,2));
xm_hat = zeros(size(rx,2),size(ry,2));
ym_hat = zeros(size(rx,2),size(ry,2));
for i = 1:size(rx,2)
    for j = 1:size(ry,2)
        xm_hat(i,:) = rx(i);
        ym_hat(:,j) = ry(j);
    end
end
xM = reshape(xm_hat,1,size(rx,2)*size(ry,2));
yM = reshape(ym_hat,1,size(rx,2)*size(ry,2));
%生成完备字典集
for i = 1:0.25*fs
    for k = 1:size(rx,2)*size(ry,2)
        Q(i,k) = (exp(-1j*4*pi*f0*((Rt+xM(k)*cosd(phi(i))+yM(k)*sind(phi(i)))/c)));
    end
end
%设置仿真时散射中心分布
a = zeros(1,size(rx,2)*size(ry,2));
a(51) = 9;
a(51*25+1) = 15;
a(end) = 13;
a(34)=10+2i;
% for i = 1:0.25*fs
%         s(i) = Q*a';
% end
%散射中心模型的回波信号（未添加噪声）
s = Q*a';
snr = 2; %修改高斯噪声的信噪比
s = awgn(s,snr,'measured'); %散射中心模型的回波信号（添加噪声）
% s1 = awgn(s1,snr,'measured'); %散射中心模型的回波信号（添加噪声）
% [ theta ] = omp( s,Q,10 ) ; %OMP算法进行散射中心估计
[theta,residual_start,residual_end] = omp2(Q, s, 10);
% plot(abs(theta));
[showindex] = find(theta~=0); %统计不为零的散射中心索引值，用于找出坐标和绘图
show_x = xM(showindex);
show_y = yM(showindex);
num = abs(theta);
show_num = reshape(num,size(rx,2),size(ry,2)); %绘制散射中心分布图
figure;
mesh(rx,ry,show_num');
xlabel('目标横坐标/m');
ylabel('目标纵坐标/m');
zlabel('散射中心系数');
grid;
legend(['SNR=',num2str(snr),'dB']);
figure;
plot(abs(theta));
find(theta>1)
