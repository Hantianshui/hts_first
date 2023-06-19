% 2023.2.13 该脚本为0.5°分辨率，采用极大似然估计每个角度的σ值，进行B矩阵的建立，
% 尾部拟合较差
% 此脚本为HMM模型运行前数据准备
%  输入：
%  O（观测序列）
%  A（状态转移矩阵）
%  B（发射矩阵）
%  PI（初始状态概率分布矩阵）
%  
%  输出：
%  估计的 A B PI矩阵，隐状态序列I 
%
clc;
close all;
clear;
 
%%%需要测试的程序内容 %%%  

C = 3e8;
theta_gate = 3;
pb = 0.94;
fre_delta = 10;
case_move = 2; % 目标运动形式选择：1，平动；2，进动
fre_num = 1; % 频点选择，1-201 9G-11G 0.01G为步长
%方差选择区间
sigmaA_line = 1:0.000001:10;
sstart = 1+7800;
eend = 500+7800;
cnt = 1;
tic 

%% 静态数据库WD数据读取
load('WD.mat','echo','ConfigWData');  % 复数数据，201*1801矩阵，横向是各个角度，纵向是各个频点

% ConfigWData.FreStart = 9;  %扫频起始频率Ghz
% ConfigWData.FreEnd = 11;   %扫频终止频率Ghz 
% ConfigWData.FreDelt = 0.01;  %扫频步进 GHz
% ConfigWData.FreNum = (ConfigWData.FreEnd - ConfigWData.FreStart)/ConfigWData.FreDelt+1;  %扫频点个数
% ConfigWData.AziStart = 0;  %方位角起始度数
% ConfigWData.AziEnd = 180;   %方位角终止度数
% ConfigWData.AziDelt = 0.1;  %方位角步进间隔
% ConfigWData.AziNum = (ConfigWData.AziEnd - ConfigWData.AziStart)/ConfigWData.AziDelt+1; %方位角个数
% NF = ConfigWData.FreNum;
% NA = ConfigWData.AziNum; 
% 
% B = (ConfigWData.FreEnd - ConfigWData.FreStart)*1e9;
% 
% figure;
% map1 = fft(echo,[],1);
% map2 = fftshift(map1,1);
% mesh(abs(map2));
% title('目标原始距离-时间像');

%% 测试数据读取
if case_move==1 % 1，平动
    % 平动情况下姿态角数据和RCS数据读取
    load('WD_pd','time','angle','echo_WD_interp');
elseif case_move==2 % 2，进动
    % 进动情况下姿态角和RCS数据读取
    load('WD_jd','time','angle','echo_WD_interp');
end

%% 读取一个频点的数据（包括观测数据和数据库数据）

angle0 = [ConfigWData.AziStart:ConfigWData.AziDelt:ConfigWData.AziEnd];%0.1°步长的0-180°
angle_len = length(angle);%angle为运动过程中角度序列，采样周期为0.1s
xx1 = 0.1*[1:angle_len];
% figure; plot(xx1,angle);
% xlabel('时间 (s)');
% ylabel('角度 (°)');
% title('真实的时间-角度');

echo_WD_interp_A  = echo_WD_interp(fre_num,:);%选定某一频点的所有采样时间RCS序列
xx2 = 0.1*[1:length(echo_WD_interp_A)];
% figure; plot(xx2,abs(echo_WD_interp_A));
% xlabel('时间 (s)');
% ylabel('幅度');
% title('测试的时间-角度');%这里应该是时间对应RCS幅度，用时间幅度去估计时间角度

echo_A = echo(fre_num,:);%数据库数据，选定频点所对应的那一行，即0-180°RCS数据
% figure; plot(0:0.1:(length(echo_A)-1)/10,abs(echo_A));
% xlabel('角度 (°)');
% ylabel('幅度');
% title('静态数据库的角度-幅度');
x = 0:0.1:180;
x1 = 0:0.01:180;
echo_A_01deg = interp1(x,echo_A,x1,'spline');
%% 将数据转化为算法所需的矩阵（找发射矩阵的分布）
angle_qua = roundn(angle,-1); %将实测的姿态角（隐状态对应于实测的RCS）进行量化，分辨率为1°

countsigma = zeros(361,1);
countsigma2 = zeros(361,1);
countmu = zeros(361,1);

%用极大似然法估计0.1°分辨率的所有均值方差
% rcsdata = abs(echo_A)/10;
% rcsdata_01deg = abs(echo_A_01deg)'/10;
rcsdata = 10*log10(abs(echo_A));
rcsdata_01deg = 10*log10(abs(echo_A_01deg)');
test_data = rcsdata_01deg((2)*50-25:(2)*50+24);
pd = zeros(1801,1);
% [paramhat,paramint]=mle(rcsdata(1:5),'distribution','norm');
[paramhat,paramint]=mle(rcsdata_01deg(1:5),'distribution','norm','Alpha',0.0001);
mu_01deg = zeros(1801,1);
sigma_01deg = zeros(1801,1);
mu_01deg(1) = paramhat(1);
sigma_01deg(1) = paramhat(2);
% sigma_1deg(1) = paramint(2,1);
for i = 2:1800
    [paramhat,paramint]=mle(rcsdata_01deg((i-1)*10-4:(i-1)*10+5),'distribution','norm','Alpha',0.0001);
    mu_01deg(i) = paramhat(1);
    sigma_01deg(i) = paramhat(2);
end
[paramhat,paramint]=mle(rcsdata_01deg(end-4:end),'distribution','norm','Alpha',0.0001);
mu_01deg(end) = paramhat(1);
sigma_01deg(end) = paramhat(2);

%————————%
%% 完全与实测姿态角无关，此条件下进行A,B,PI矩阵的假设
%隐状态数目构造初始状态转移矩阵A
N_I = 1801;
%选择观测序列
% O = abs(echo_WD_interp_A(1+mark:long+mark));
% O = abs(echo_WD_interp_A(1:end));
O = 10*log10(abs(echo_WD_interp_A(sstart:eend)));
%用来做对比的角度值（已量化）
% angle_com = angle_qua(1+mark:long+mark);
angle_com = angle_qua(sstart:eend);
%得到所有的不重复含误差的观测状态,且从小到大排序
Oobsv_unq = unique(O); 

%将观测序列由数据转化为索引值
for i = 1:length(O)
    O_num(i) = find(Oobsv_unq==O(i));
end

%初始化A、B矩阵
A = zeros(N_I);
%计算隐状态对应的数据库RCS与观测RCS的误差，代入概率密度函数求概率，得到发射矩阵
B = zeros(N_I,length(Oobsv_unq));
%发射矩阵每一行满足正态分布，其中参数已经通过极大似然估计得到
%构造B矩阵
for i = 1:N_I
        sigma = sigma_01deg(i);
        mu_O = mu_01deg(i);
    for j = 1:length(Oobsv_unq)
        B(i,j) =  normpdf(Oobsv_unq(j),mu_O,sigma);%使用高斯分布构造
    end
end

%将B矩阵按每一行变成概率和为1
B = B./repmat(sum(B,2),1,size(B,2));
B(isnan(B)) = 0;
% 初始状态概率分布矩阵(平均分布)
PI = ones(N_I,1);
% PI(1202)=10000;
PI = PI./repmat(sum(PI),N_I,1);

%循环搜索找到最合适的A矩阵方差，然后退出搜索
%给定最大匹配概率值和A矩阵方差初值
max_true = 0;
% sigmaA = 2;
max_dist = 0;
% sigma_delta = 0.02;

%A矩阵参数迭代
while(1)
sigmaA = sigmaA_line(min(find(normcdf(theta_gate,0,sigmaA_line)<=pb)));
%构造A矩阵，使用高斯分布，并迭代减小方差值
for i = 1:N_I
    for j = 1:N_I
        A(i,j) =  normpdf(j,i,sigmaA);
    end
end
%将A矩阵按每一行变成概率和为1
A = A./repmat(sum(A,2),1,size(A,2));
A(isnan(A)) = 0;



%% 进行隐状态序列的估计
% 输入要求，
% O_num表示观测状态的索引 为N_I*1 列向量
% A表示状态转移矩阵 为N*N维 方阵
% B表示发射矩阵 为N*M维 行数表示隐状态数目，列数表示观测状态数目
% Pi表示初始隐状态概率分布矩阵 为N*1维 列向量
A_GUESS = A;
B_GUESS = B;
PI_GUESS = PI;
%viterbi算法
% [guessTR,guessE,guessPi,logliks] = hmmtrain_system(O_num,A_GUESS,B_GUESS,PI_GUESS,'TOLERANCE',2.2e-4);
%近似算法
% [~,logPseq,fs,bs,scale] = hmmdecode_system(O_num,guessTR,guessE,guessPi);
% f = fs.*repmat(cumprod(scale),size(fs,1),1);
% [~,index] = max(fs.*bs,[],1);
% indexangle = ((index-1)*0.5);
% hold on;
% %估计值
% plot(indexangle);
% [likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, guessTR,guessE,guessPi);
[likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, A_GUESS,B_GUESS,PI_GUESS);
% [likelystates,logp,pTR,vall] = hmmviterbi_system2(O_num, guessTR,guessE,guessPi);
% filter算法
% for i = 2:length(O_num)-100+1
%     filterseq = O_num(i:i+100-1);
%     [guessTR,guessE,guessPi,logliks] = hmmtrain_system(filterseq,A_hat,B_hat,Pi,'TOLERANCE',1e-3);
% %     [pStates,pSeq, fs, bs, s_ori] = hmmdecode_system(O_num,guessTR,guessE,guessPi);
% %     filterseq = exp(log(fs)+log(bs));
% %     [~,I_filter] = max(filterseq);
% end


%% 绘制结果
%隐状态[1-1801]与真实的姿态角[0-180]相差-1/10
likelyangle = ((likelystates-1)/10);
prob = sum(angle_com' == likelyangle)/length(likelyangle);
% a = DiscreteFrechetDist(likelyangle,angle_com')
%每次保留最大值，当出现峰值时，退出迭代
% if prob > max_true
%     max_true = prob;
%     fprintf('目前最佳概率%8.4f\n',prob);
%     sigmaA = sigmaA - sigma_delta;
% else
%     if a > max_dist
%         max_dist = a;
%     else
%         break;
%     end
%     break;
% end


%用作对比的真值
prob = prob*100;
% figure;
leg_str{cnt} = [num2str(theta_gate),'°','内概率为',num2str(pb),'估计值'];  
lastwork(:,cnt) = likelyangle;
cnt = cnt + 1;


%估计值
% plot(likelyangle);
frq = (fre_num-1)*0.01+9;
% hold on;
% title(['频率为',num2str(frq),'GHz时,方差为：',num2str(sigmaA),',准确概率为：',num2str(prob),'%']);
% figure;
% % title('pdf图');
% Delta = (angle_com' - likelyangle);
% cdfplot(Delta);
toc  
disp(['频率为',num2str(frq),'GHz时，运行时间为：',num2str(toc)]);

% if fre_num < 192
%     fre_num = fre_num + fre_delta;
% else
%     break;
% end
if pb < 0.98
    pb = pb + 0.01;
    
else
%     plot(angle_com);
    leg_str{cnt} = '真值';
%     legend(leg_str);
    break;
end


end

figure;
plot(lastwork);
grid on;
hold on;
plot(angle_com);
legend(leg_str);
title([num2str(frq),'GHz时,',num2str(sstart),'-',num2str(eend),'序列估计结果']);

