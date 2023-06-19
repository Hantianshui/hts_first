
% 2023.2.13 �ýű�Ϊ0.5��ֱ��ʣ����ü�����Ȼ����ÿ���ǶȵĦ�ֵ������B����Ľ�����
% β����Ͻϲ�
% �˽ű�ΪHMMģ������ǰ����׼��
%  ���룺
%  O���۲����У�
%  A��״̬ת�ƾ���
%  B���������
%  PI����ʼ״̬���ʷֲ�����
%  
%  �����
%  ���Ƶ� A B PI������״̬����I 
%
clc;
close all;
clear;
 
%%%��Ҫ���Եĳ������� %%%  

C = 3e8;
sigmaA_line = 1:0.001:2;
% sigmaA = sigmaA_line(min(find(normpdf(71,1,sigmaA_line)>0)));
theta_gate = 3;
pb = 0.94;
sigmaA = sigmaA_line(min(find(normcdf(theta_gate,0,sigmaA_line)<=pb)));
fre_delta = 10;
case_move = 2; % Ŀ���˶���ʽѡ��1��ƽ����2������
fre_num = 1; % Ƶ��ѡ��1-201 9G-11G 0.01GΪ����
while(1)
tic 
%% ��̬���ݿ�WD���ݶ�ȡ
load('WD.mat','echo','ConfigWData');  % �������ݣ�201*1801���󣬺����Ǹ����Ƕȣ������Ǹ���Ƶ��

% ConfigWData.FreStart = 9;  %ɨƵ��ʼƵ��Ghz
% ConfigWData.FreEnd = 11;   %ɨƵ��ֹƵ��Ghz 
% ConfigWData.FreDelt = 0.01;  %ɨƵ���� GHz
% ConfigWData.FreNum = (ConfigWData.FreEnd - ConfigWData.FreStart)/ConfigWData.FreDelt+1;  %ɨƵ�����
% ConfigWData.AziStart = 0;  %��λ����ʼ����
% ConfigWData.AziEnd = 180;   %��λ����ֹ����
% ConfigWData.AziDelt = 0.1;  %��λ�ǲ������
% ConfigWData.AziNum = (ConfigWData.AziEnd - ConfigWData.AziStart)/ConfigWData.AziDelt+1; %��λ�Ǹ���
% NF = ConfigWData.FreNum;
% NA = ConfigWData.AziNum; 
% 
% B = (ConfigWData.FreEnd - ConfigWData.FreStart)*1e9;
% 
% figure;
% map1 = fft(echo,[],1);
% map2 = fftshift(map1,1);
% mesh(abs(map2));
% title('Ŀ��ԭʼ����-ʱ����');

%% �������ݶ�ȡ
if case_move==1 % 1��ƽ��
    % ƽ���������̬�����ݺ�RCS���ݶ�ȡ
    load('WD_pd','time','angle','echo_WD_interp');
elseif case_move==2 % 2������
    % �����������̬�Ǻ�RCS���ݶ�ȡ
    load('WD_jd','time','angle','echo_WD_interp');
end

%% ��ȡһ��Ƶ������ݣ������۲����ݺ����ݿ����ݣ�

angle0 = [ConfigWData.AziStart:ConfigWData.AziDelt:ConfigWData.AziEnd];%0.1�㲽����0-180��
angle_len = length(angle);%angleΪ�˶������нǶ����У���������Ϊ0.1s
xx1 = 0.1*[1:angle_len];
% figure; plot(xx1,angle);
% xlabel('ʱ�� (s)');
% ylabel('�Ƕ� (��)');
% title('��ʵ��ʱ��-�Ƕ�');

echo_WD_interp_A  = echo_WD_interp(fre_num,:);%ѡ��ĳһƵ������в���ʱ��RCS����
xx2 = 0.1*[1:length(echo_WD_interp_A)];
% figure; plot(xx2,abs(echo_WD_interp_A));
% xlabel('ʱ�� (s)');
% ylabel('����');
% title('���Ե�ʱ��-�Ƕ�');%����Ӧ����ʱ���ӦRCS���ȣ���ʱ�����ȥ����ʱ��Ƕ�

echo_A = echo(fre_num,:);%���ݿ����ݣ�ѡ��Ƶ������Ӧ����һ�У���0-180��RCS����
% figure; plot(0:0.1:(length(echo_A)-1)/10,abs(echo_A));
% xlabel('�Ƕ� (��)');
% ylabel('����');
% title('��̬���ݿ�ĽǶ�-����');
x = 0:0.1:180;
x1 = 0:0.01:180;
echo_A_01deg = interp1(x,echo_A,x1,'spline');
%% ������ת��Ϊ�㷨����ľ����ҷ������ķֲ���
angle_qua = roundn(angle,-1); %��ʵ�����̬�ǣ���״̬��Ӧ��ʵ���RCS�������������ֱ���Ϊ1��

% deltasigma = zeros(361,1);
% deltamu = zeros(361,1);
countsigma = zeros(361,1);
countsigma2 = zeros(361,1);
countmu = zeros(361,1);
% for i = 1:361
%     ind = find(angle_qua>i-0.25&angle_qua<=i+0.25);
% %     deltadata = delta(ind);
%     countrcs = abs(echo_WD_interp_A(ind));
%     countdelta = max(countrcs)-min(countrcs);
%     if(isempty(countrcs)==0)
% %         [paramhat,~] = mle(deltadata,'distribution','norm');
% %         deltamu(i+1) = paramhat(1);
% %         deltasigma(i+1) = paramhat(2);
%         [paramhat,~] = mle(countrcs,'distribution','norm');
%         countmu(i+1) = paramhat(1);
%         countsigma(i+1) = paramhat(2);
%         countsigma2(i+1) = countdelta;
% %         deltadata = [];
%         countrcs = [];
%         countdelta = [];
%         ind = [];
%     end
% end


%�ü�����Ȼ������0.1��ֱ��ʵ����о�ֵ����
% rcsdata = abs(echo_A)/10;
% rcsdata_01deg = abs(echo_A_01deg)'/10;
rcsdata = 10*log10(abs(echo_A));
rcsdata_01deg = 10*log10(abs(echo_A_01deg)');
test_data = rcsdata_01deg((2)*50-25:(2)*50+24);
pd = zeros(1801,1);
% [paramhat,paramint]=mle(rcsdata(1:5),'distribution','norm');
[paramhat,paramint]=mle(rcsdata_01deg(1:5),'distribution','norm','Alpha',0.0001);
% [phat,pci] = mle(rcsdata_01deg(1:5),'pdf',@(x,s,sigma) pdf('rician',x,s,5),'start',10);
mu_01deg = zeros(1801,1);
sigma_01deg = zeros(1801,1);
mu_01deg(1) = paramhat(1);
sigma_01deg(1) = paramhat(2);
% sigma_1deg(1) = paramint(2,1);
for i = 2:1800
%     [paramhat,paramint]=mle(rcsdata((i-1)*10-4:(i-1)*10+5),'distribution','norm');
    [paramhat,paramint]=mle(rcsdata_01deg((i-1)*10-4:(i-1)*10+5),'distribution','norm','Alpha',0.0001);
%     if i<60*2
%         sigma_05deg(i) = paramhat(2);
%         mu_05deg(i) = paramhat(1);
%     else
        mu_01deg(i) = paramhat(1);
        sigma_01deg(i) = paramhat(2);
%     end
%     sigma_1deg(i) = paramint(2,1);
end
% [paramhat,paramint]=mle(rcsdata(end-4:end),'distribution','norm');
[paramhat,paramint]=mle(rcsdata_01deg(end-4:end),'distribution','norm','Alpha',0.0001);
mu_01deg(end) = paramhat(1);
sigma_01deg(end) = paramhat(2);

%����������������%

%% ��ȫ��ʵ����̬���޹أ��������½���A,B,PI����ļ���
%��״̬��Ŀ�����ʼ״̬ת�ƾ���A
N_I = 1801;
%ѡ��۲�����
% O = abs(echo_WD_interp_A(1+mark:long+mark));
% O = abs(echo_WD_interp_A(1:end));
O = 10*log10(abs(echo_WD_interp_A(1:end)));
%�������ԱȵĽǶ�ֵ����������
% angle_com = angle_qua(1+mark:long+mark);
angle_com = angle_qua(1:end);
%�õ����еĲ��ظ������Ĺ۲�״̬,�Ҵ�С��������
Oobsv_unq = unique(O); 

%���۲�����������ת��Ϊ����ֵ
for i = 1:length(O)
    O_num(i) = find(Oobsv_unq==O(i));
end

%��ʼ��A��B����
A = zeros(N_I);
%������״̬��Ӧ�����ݿ�RCS��۲�RCS������������ܶȺ�������ʣ��õ��������
B = zeros(N_I,length(Oobsv_unq));
%�������ÿһ��������̬�ֲ������в����Ѿ�ͨ��������Ȼ���Ƶõ�
%����B����
for i = 1:N_I
        sigma = sigma_01deg(i);
        mu_O = mu_01deg(i);
%     if i == 1 
%         pd = fitdist(rcsdata_01deg(1:5),'Rician');
%     elseif i >= 2 && i<=1800
%         pd = fitdist(rcsdata_01deg((i-1)*10-4:(i-1)*10+5),'Rician');
%     else
%         pd = fitdist(rcsdata_01deg(end-4:end),'Rician');
%     end
    for j = 1:length(Oobsv_unq)
%         if Oobsv_unq(j)>= min_mle(i)&&Oobsv_unq(j)<= max_mle(i)
        B(i,j) =  normpdf(Oobsv_unq(j),mu_O,sigma);%ʹ�ø�˹�ֲ�����
    end
end

%��B����ÿһ�б�ɸ��ʺ�Ϊ1
B = B./repmat(sum(B,2),1,size(B,2));
B(isnan(B)) = 0;
% ��ʼ״̬���ʷֲ�����(ƽ���ֲ�)
PI = ones(N_I,1);
PI = PI./repmat(sum(PI),N_I,1);

%ѭ�������ҵ�����ʵ�A���󷽲Ȼ���˳�����
%�������ƥ�����ֵ��A���󷽲��ֵ
max_true = 0;
% sigmaA = 2;
max_dist = 0;
% sigma_delta = 0.02;

%����A����ʹ�ø�˹�ֲ�����������С����ֵ
for i = 1:N_I
    for j = 1:N_I
%         if i>=j
%             A(i,j) = r^(abs(i-j));
%         else
%             A(i,j) = r^(abs(i-j));
%         end
        A(i,j) =  normpdf(j,i,sigmaA);

    end
end
%��A����ÿһ�б�ɸ��ʺ�Ϊ1
% A(find(A<r^(60))) = 0;
A = A./repmat(sum(A,2),1,size(A,2));
A(isnan(A)) = 0;



%% ������״̬���еĹ���
% ����Ҫ��
% O_num��ʾ�۲�״̬������ ΪN_I*1 ������
% A��ʾ״̬ת�ƾ��� ΪN*Nά ����
% B��ʾ������� ΪN*Mά ������ʾ��״̬��Ŀ��������ʾ�۲�״̬��Ŀ
% Pi��ʾ��ʼ��״̬���ʷֲ����� ΪN*1ά ������
A_GUESS = A;
B_GUESS = B;
PI_GUESS = PI;
%viterbi�㷨
% [guessTR,guessE,guessPi,logliks] = hmmtrain_system(O_num,A_GUESS,B_GUESS,PI_GUESS,'TOLERANCE',2.2e-4);
%�����㷨
% [~,logPseq,fs,bs,scale] = hmmdecode_system(O_num,guessTR,guessE,guessPi);
% f = fs.*repmat(cumprod(scale),size(fs,1),1);
% [~,index] = max(fs.*bs,[],1);
% indexangle = ((index-1)*0.5);
% hold on;
% %����ֵ
% plot(indexangle);
% [likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, guessTR,guessE,guessPi);
[likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, A_GUESS,B_GUESS,PI_GUESS);
% [likelystates,logp,pTR,vall] = hmmviterbi_system2(O_num, guessTR,guessE,guessPi);
% filter�㷨
% for i = 2:length(O_num)-100+1
%     filterseq = O_num(i:i+100-1);
%     [guessTR,guessE,guessPi,logliks] = hmmtrain_system(filterseq,A_hat,B_hat,Pi,'TOLERANCE',1e-3);
% %     [pStates,pSeq, fs, bs, s_ori] = hmmdecode_system(O_num,guessTR,guessE,guessPi);
% %     filterseq = exp(log(fs)+log(bs));
% %     [~,I_filter] = max(filterseq);
% end


%% ���ƽ��
%��״̬[1-1801]����ʵ����̬��[0-180]���-1/10
likelyangle = ((likelystates-1)/10);
prob = sum(angle_com' == likelyangle)/length(likelyangle);
a = DiscreteFrechetDist(likelyangle,angle_com')
%ÿ�α������ֵ�������ַ�ֵʱ���˳�����
% if prob > max_true
%     max_true = prob;
%     fprintf('Ŀǰ��Ѹ���%8.4f\n',prob);
%     sigmaA = sigmaA - sigma_delta;
% else
%     if a > max_dist
%         max_dist = a;
%     else
%         break;
%     end
%     break;
% end


%�����Աȵ���ֵ
prob = prob*100;
figure;
plot(angle_com);
hold on;
%����ֵ
plot(likelyangle);
legend('��ֵ','����ֵ');
frq = (fre_num-1)*0.01+9;
title(['Ƶ��Ϊ',num2str(frq),'GHzʱ���ƽ��']);
xlabel('���к�');
ylabel('�Ƕ�(��)');
% figure;
% % title('pdfͼ');
% Delta = (angle_com' - likelyangle);
% cdfplot(Delta);
toc  
disp(['Ƶ��Ϊ',num2str(frq),'GHzʱ������ʱ��Ϊ��',num2str(toc)]);

if fre_num < 192
    fre_num = fre_num + fre_delta;
else
    break;
end

end
