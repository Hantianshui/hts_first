%%��������������������������������������������%%
% 2023.2.22 �˽ű�Ϊ1��ֱ��ʲ��ü�����Ȼ���Ʒ��չʾ�����뼫����Ȼ���ƽ��������Ҫ����������ݿ�����Ӷϵ�
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

C = 3e8;

case_move = 2; % Ŀ���˶���ʽѡ��1��ƽ����2������
fre_num = 101; % Ƶ��ѡ��1-201 9G-11G 0.01GΪ����
%% ��̬���ݿ�WD���ݶ�ȡ
load('WD.mat','echo','ConfigWData');  % �������ݣ�201*1801���󣬺����Ǹ����Ƕȣ������Ǹ���Ƶ��

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
%��ֵ�������
x = 0:0.1:180;
x1 = 0:0.01:180;
echo_A_01deg = interp1(x,echo_A,x1,'pchip');
%% ������ת��Ϊ�㷨����ľ����ҷ������ķֲ���
angle_qua = roundn(angle,0); %��ʵ�����̬�ǣ���״̬��Ӧ��ʵ���RCS�������������ֱ���Ϊ1��
% figure; plot(xx1,angle);
%�Ƕ������ݿ��ж�Ӧ��index
angle_index = angle_qua*10+1;
%ת��Ϊ���ݿ��е�ֵ�����������������ֲ�
echo_WD = echo_A(angle_index);
%�������������������������ֲ�
delta = -abs(echo_WD(:))+abs(echo_WD_interp_A(:)); 

% rcsdata = abs(echo_A);%��������ģ
count_gate = zeros(181,1);
rcsdata = 10*log10(abs(echo_A))';
% rcsdata_01deg = power(10,abs(echo_A_01deg)'/10);

%�ü�����Ȼ������1��ֱ��ʵ����о�ֵ����
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



%% %����������������%
%ͳ�Ʒ����ó�A����
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


%% ��ȫ��ʵ����̬���޹أ��������½���A,B,PI����ļ���
%��״̬��Ŀ�����ʼ״̬ת�ƾ���A
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
%��A����ÿһ�б�ɸ��ʺ�Ϊ1
A = A./repmat(sum(A,2),1,size(A,2));
% A(find(A< r^(6))) = 0;
%NaNֵ����
A(isnan(A)) = 0;

STATE = 1;
%ѡ��۲�����
% O = abs(echo_WD_interp_A(STATE:end));
O = 10*log10(abs(echo_WD_interp_A(STATE:end)));
%�������ԱȵĽǶ�ֵ����������
angle_com = angle_qua(STATE:end);

%�۲�״̬181�� ����ֵ-1Ϊʵ�ʽǶ�
RCSDATA = echo_A(1:10:end);
% echo_A_1deg = abs(echo_A(1:10:end));
echo_A_1deg = 10*log10(abs(echo_A(1:10:end)));
%�õ����е����ݿ��������Ĺ۲�״̬,�Ҵ�С��������
Odata_unq = unique(echo_A_1deg); 
%�õ����еĲ��ظ������Ĺ۲�״̬,�Ҵ�С��������
Oobs_unq = unique(O); 
%���۲�����������ת��Ϊ����ֵ
for i = 1:length(O)
    O_num(i) = find(Oobs_unq==O(i));
end


%% ���޻���
%ͳ��0-180��ÿ���Ƕȵ������Сֵ
O_index = zeros(N_I,length(O_num));
cnt = ones(N_I,1);

for j = 1:length(O_num)
    O_index(angle_com(j)+1,cnt(angle_com(j)+1)) = O_num(j);
    cnt(angle_com(j)+1) = cnt(angle_com(j)+1)+1;
end


%ͨ��ͳ�Ƶõ�ÿ���Ƕȵ������С����ֵ
max_deg = zeros(181,1);
min_deg = zeros(181,1);
%������Ȼ���Ƶõ�ÿ���Ƕȵ������С����ֵ
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
% tilte('��ϵķֲ�Ϊ��normal');

%��������ͼ
% figure;
% plot(echo_A_1deg(10:70),'b-'); %���ݿ�����
% hold on;
% % plot(mu_1deg,'r-');
% plot(max_deg(10:70),'r.'); %ͳ������ÿ���Ƕ����ֵ
% plot(min_deg(10:70),'black.'); %ͳ������ÿ���Ƕ���Сֵ
% title('10-70������ֵ�ֲ�ͼ��');
% plot(min_mle,'r-'); %������Ȼ��������ÿ���Ƕ���Сֵ
% plot(max_mle,'black-');%������Ȼ��������ÿ���Ƕ����ֵ
% max_mle(74) = 0;
% max_mle(75) = 0;
% max_mle(92) = 0;
% min_mle(74) = 0;
% min_mle(75) = 0;
% min_mle(92) = 0;

%% ����B����
%������״̬��Ӧ�����ݿ�RCS��۲�RCS������������ܶȺ�������ʣ��õ��������
B = zeros(N_I,length(Oobs_unq));
G = zeros(N_I,length(Oobs_unq));


%���ɹ۲����B����ʱʹ��ͳ�Ƶ�����ֵ��
for i = 1:N_I
    mu = mu_1deg(i);
    sigma = sigma_1deg(i);
    %�����ܶȴ�����ʣ�����ÿһ����״̬���������й۲�״̬
    for j = 1:length(Oobs_unq)
%         if Oobs_unq(j)>= min_mle(i)&&Oobs_unq(j)<= max_mle(i)
            G(i,j) =  normpdf(Oobs_unq(j),mu,sigma);%ʹ�ø�˹�ֲ������
%             G(i,j) = normpdf(mu,mu,sigma);
%           end
%             G(i,j) = 10;
%         end
    end
end

%����ʹ�õ�B����
B = G;

%��B����ÿһ�б�ɸ��ʺ�Ϊ1
B = B./repmat(sum(B,2),1,size(B,2));
B(isnan(B)) = 0;

% ��ʼ״̬���ʷֲ�����(ƽ���ֲ�)
PI = ones(N_I,1);
% PI = zeros(N_I,1);
% PI(134) = 1;
PI = PI./repmat(sum(PI),N_I,1);


%% ������״̬���еĹ���
% ����Ҫ��
% O_num��ʾ�۲�״̬������ ΪN_I*1 ������
% A��ʾ״̬ת�ƾ��� ΪN*Nά ����
% B��ʾ������� ΪN*Mά ������ʾ��״̬��Ŀ��������ʾ�۲�״̬��Ŀ
% Pi��ʾ��ʼ��״̬���ʷֲ����� ΪN*1ά ������

A_GUESS = A;
B_GUESS = B;
PI_GUESS = PI;
% load('abc.mat');
% A_GUESS = guessTR;
% B_GUESS = guessE;
% PI_GUESS = guessPi;
%viterbi�㷨
[guessTR,guessE,guessPi,logliks] = hmmtrain_system2(O_num,O,mu_1deg,A_GUESS,B_GUESS,PI_GUESS,'maxiterations',1);
% [likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, A_GUESS,B_GUESS,PI_GUESS);
[likelystates,logp,pTR,vall] = hmmviterbi_system(O_num, guessTR,guessE,guessPi);
% % �����㷨
% hmmtrain_system
% [~,logPseq,fs,bs,scale] = hmmdecode_system(O_num, A_GUESS,B_GUESS,PI_GUESS);
% % f = fs.*repmat(cumprod(scale),size(fs,1),1);
% [~,index] = max(fs.*bs,[],1);
% indexangle = ((index-1));
% save('abc.mat','guessTR','guessE','guessPi');

%% ���ƽ��
%��״̬[1-181]����ʵ����̬��[0-180]���1
likelyangle = (likelystates-1);
% likelyangle = round(smooth(likelyangle,10))';
a = sum(angle_com' == likelyangle)/length(likelyangle)*100;
%�����Աȵ���ֵ
figure;
plot(angle_com);
hold on;
%����ֵ
plot(likelyangle);
legend('��ֵ','����ֵ');
title(['׼ȷ����Ϊ��',num2str(a),'%']);
figure;
% title('pdfͼ');
Delta = angle_com' - likelyangle;
cdfplot(Delta);
DiscreteFrechetDist(likelyangle,angle_com')