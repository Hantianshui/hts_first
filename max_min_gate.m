%%——————————————————————%%
% 2023.2.22 此脚本为1°分辨率采用极大似然估计方差，展示门限与极大似然估计结果，不需要最终拟合数据可以添加断点
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

C = 3e8;

case_move = 2; % 目标运动形式选择：1，平动；2，进动
fre_num = 101; % 频点选择，1-201 9G-11G 0.01G为步长
%% 静态数据库WD数据读取
load('WD.mat','echo','ConfigWData');  % 复数数据，201*1801矩阵，横向是各个角度，纵向是各个频点

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
%插值后的数据
x = 0:0.1:180;
x1 = 0:0.01:180;
echo_A_01deg = interp1(x,echo_A,x1,'pchip');
%% 将数据转化为算法所需的矩阵（找发射矩阵的分布）
angle_qua = roundn(angle,0); %将实测的姿态角（隐状态对应于实测的RCS）进行量化，分辨率为1°
% figure; plot(xx1,angle);
%角度在数据库中对应的index
angle_index = angle_qua*10+1;
%转化为数据库中的值，用来分析测量误差分布
echo_WD = echo_A(angle_index);
%整体数据误差，可以用来分析误差分布
delta = -abs(echo_WD(:))+abs(echo_WD_interp_A(:)); 

% rcsdata = abs(echo_A);%将复数求模
count_gate = zeros(181,1);
rcsdata = 10*log10(abs(echo_A))';
% rcsdata_01deg = power(10,abs(echo_A_01deg)'/10);

%用极大似然法估计1°分辨率的所有均值方差
mu_1deg = zeros(181,1);
sigma_1deg = zeros(181,1);
[paramhat,paramint]=mle(rcsdata(1:5),'distribution','Normal','alpha',0.0001);
% count_gate(1) = max(rcsdata(1:5))-min(rcsdata(1:5));
% [paramhat,paramint]=mle(rcsdata(1:5),'distribution','norm');
mu_1deg(1) = paramhat(1);
% sigma_1deg(1) = (paramint(2,1)+paramint(2,2))/2;
sigma_1deg(1) = paramhat(2);
for i = 2:180
    [paramhat,paramint]=mle(rcsdata((i-1)*10-4:(i-1)*10+5),'distribution','Normal','alpha',0.0001);
    mu_1deg(i) = paramhat(1);
%     sigma_1deg(i) = (paramint(2,1)+paramint(2,2))/2;
%     count_gate(i) = max(rcsdata((i-1)*10-4:(i-1)*10+5))-min(rcsdata((i-1)*10-4:(i-1)*10+5));
    sigma_1deg(i) = paramhat(2);
end
[paramhat,paramint]=mle(rcsdata(end-4:end),'distribution','Normal','alpha',0.0001);
mu_1deg(end) = paramhat(1);
% sigma_1deg(end) = (paramint(2,1)+paramint(2,2))/2;
sigma_1deg(end) = paramhat(2);
% count_gate(end) = max(rcsdata(end-4:end))-min(rcsdata(end-4:end));

% cnt_gate = sort(count_gate(1:80));
% sigma_1 = median(cnt_gate)/6;
% cnt_gate = sort(count_gate(81:90));
% sigma_2 = median(cnt_gate)/6;
% cnt_gate = sort(count_gate(91:end));
% sigma_3 = median(cnt_gate)/6;



%% %————————%
%统计方法得出A矩阵
C = zeros(length(angle_qua)-1,2);
for i = 1:size(C,1)
    temp = [angle_qua(i),angle_qua(i+1)];
    C(i,:) = temp;
end
[b,m,n] = unique(C,'rows');
c = tabulate(n);
cnt = m(c(:,1));
s = sortrows([C(m,:),c(:,2:end)],3);
A_count = zeros(181,181);
I_uni = 0:1:180;
for cnt = 1:size(s,1)
    i = find(I_uni==s(cnt,1));
    j = find(I_uni==s(cnt,2));
    A_count(i,j) = s(cnt,3);
end
A_count = A_count./repmat(sum(A_count,2),1,size(A_count,2));
A_count(isnan(A_count)) = 0;


%% 完全与实测姿态角无关，此条件下进行A,B,PI矩阵的假设
%隐状态数目构造初始状态转移矩阵A
N_I = 181;
A = zeros(N_I);
r = 0.03;
for i = 1:N_I
    for j = 1:N_I
%         if i>=j
%             A(i,j) = r^(abs(i-j));
%         else
%             A(i,j) = r^(abs(i-j));
%         end
        A(i,j) =  normpdf(j,i,0.4);
    end
end
%将A矩阵按每一行变成概率和为1
A = A./repmat(sum(A,2),1,size(A,2));
% A(find(A< r^(6))) = 0;
%NaN值置零
A(isnan(A)) = 0;

STATE = 1;
%选择观测序列
% O = abs(echo_WD_interp_A(STATE:end));
O = 10*log10(abs(echo_WD_interp_A(STATE:end)));
%用来做对比的角度值（已量化）
angle_com = angle_qua(STATE:end);

%观测状态181个 索引值-1为实际角度
RCSDATA = echo_A(1:10:end);
% echo_A_1deg = abs(echo_A(1:10:end));
echo_A_1deg = 10*log10(abs(echo_A(1:10:end)));
%得到所有的数据库中无误差的观测状态,且从小到大排序
Odata_unq = unique(echo_A_1deg); 
%得到所有的不重复含误差的观测状态,且从小到大排序
Oobs_unq = unique(O); 
%将观测序列由数据转化为索引值
for i = 1:length(O)
    O_num(i) = find(Oobs_unq==O(i));
end


%% 门限绘制
%统计0-180°每个角度的最大最小值
O_index = zeros(N_I,length(O_num));
cnt = ones(N_I,1);

for j = 1:length(O_num)
    O_index(angle_com(j)+1,cnt(angle_com(j)+1)) = O_num(j);
    cnt(angle_com(j)+1) = cnt(angle_com(j)+1)+1;
end


%通过统计得到每个角度的最大最小门限值
max_deg = zeros(181,1);
min_deg = zeros(181,1);
%极大似然估计得到每个角度的最大最小门限值
max_mle = zeros(181,1);
min_mle = zeros(181,1);
M = 360;
RCS_deg = zeros(181,M);
m = 1;
for i = 1:181
    if(O_index(i,1)~=0)
    max_deg(i) = Oobs_unq(max(O_index(i,find(O_index(i,:)>0))));
    min_deg(i) = Oobs_unq(min(O_index(i,find(O_index(i,:)>0))));
    a = Oobs_unq((O_index(i,find(O_index(i,:)>0))));
    a=[a , zeros(1,M-length(a))];
    RCS_deg(i,:) = a;
    end
%     if i>83&&i<87
%         max_mle(i) = mu_1deg(i)+2;
%         min_mle(i) = mu_1deg(i)-2;
%     else
    K = 1.5;
    max_mle(i) = mu_1deg(i) + K*sigma_1deg(i);
    min_mle(i) = mu_1deg(i) - K*sigma_1deg(i);
%     end
end
% for i = 16:1:134
%     K = i;
%     RCS = RCS_deg(K,(find(RCS_deg(K,:))))-echo_A_1deg(i);
%     histfit(RCS,5);
% end
% RCS = cell2mat(RCS_deg);
% max(RCS_deg);

% 
% RCS = RCS_deg(130,(find(RCS_deg(130,:))));
% figure;
% h = histfit(RCS,18,'weibull');
% tilte('拟合的分布为：normal');

%绘制门限图
% figure;
% plot(echo_A_1deg(10:70),'b-'); %数据库数据
% hold on;
% % plot(mu_1deg,'r-');
% plot(max_deg(10:70),'r.'); %统计所得每个角度最大值
% plot(min_deg(10:70),'black.'); %统计所得每个角度最小值
% title('10-70°门限值分布图像');
% plot(min_mle,'r-'); %极大似然估计所得每个角度最小值
% plot(max_mle,'black-');%极大似然估计所得每个角度最大值
% max_mle(74) = 0;
% max_mle(75) = 0;
% max_mle(92) = 0;
% min_mle(74) = 0;
% min_mle(75) = 0;
% min_mle(92) = 0;

%% 生成B矩阵
%计算隐状态对应的数据库RCS与观测RCS的误差，代入概率密度函数求概率，得到发射矩阵
B = zeros(N_I,length(Oobs_unq));
G = zeros(N_I,length(Oobs_unq));


%生成观测矩阵B（此时使用统计的门限值）
for i = 1:N_I
    mu = mu_1deg(i);
    sigma = sigma_1deg(i);
    %概率密度代替概率，对于每一个隐状态，遍历所有观测状态
    for j = 1:length(Oobs_unq)
%         if Oobs_unq(j)>= min_mle(i)&&Oobs_unq(j)<= max_mle(i)
            G(i,j) =  normpdf(Oobs_unq(j),mu,sigma);%使用高斯分布求概率
%             G(i,j) = normpdf(mu,mu,sigma);
%           end
%             G(i,j) = 10;
%         end
    end
end

%代入使用的B矩阵
B = G;

%将B矩阵按每一行变成概率和为1
B = B./repmat(sum(B,2),1,size(B,2));
B(isnan(B)) = 0;

% 初始状态概率分布矩阵(平均分布)
PI = ones(N_I,1);
% PI = zeros(N_I,1);
% PI(134) = 1;
PI = PI./repmat(sum(PI),N_I,1);


%% 进行隐状态序列的估计
% 输入要求，
% O_num表示观测状态的索引 为N_I*1 列向量
% A表示状态转移矩阵 为N*N维 方阵
% B表示发射矩阵 为N*M维 行数表示隐状态数目，列数表示观测状态数目
% Pi表示初始隐状态概率分布矩阵 为N*1维 列向量

A_GUESS = A;
B_GUESS = B;
PI_GUESS = PI;
% load('abc.mat');
% A_GUESS = guessTR;
% B_GUESS = guessE;
% PI_GUESS = guessPi;
%viterbi算法
[guessTR,guessE,guessPi,logliks] = hmmtrain_system2(O_num,O,mu_1deg,A_GUESS,B_GUESS,PI_GUESS,'maxiterations',1);
% [likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, A_GUESS,B_GUESS,PI_GUESS);
[likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, guessTR,guessE,guessPi);
% % 近似算法
% hmmtrain_system
% [~,logPseq,fs,bs,scale] = hmmdecode_system(O_num, A_GUESS,B_GUESS,PI_GUESS);
% % f = fs.*repmat(cumprod(scale),size(fs,1),1);
% [~,index] = max(fs.*bs,[],1);
% indexangle = ((index-1));
% save('abc.mat','guessTR','guessE','guessPi');

%% 绘制结果
%隐状态[1-181]与真实的姿态角[0-180]相差1
likelyangle = (likelystates-1);
% likelyangle = round(smooth(likelyangle,10))';
a = sum(angle_com' == likelyangle)/length(likelyangle)*100;
%用作对比的真值
figure;
plot(angle_com);
hold on;
%估计值
plot(likelyangle);
legend('真值','估计值');
title(['准确概率为：',num2str(a),'%']);
figure;
% title('pdf图');
Delta = angle_com' - likelyangle;
cdfplot(Delta);
DiscreteFrechetDist(likelyangle,angle_com')