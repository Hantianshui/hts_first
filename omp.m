function [x_hat,residual_start,residual_end] = omp(D, y, K)
    % 输入参数：
    % D：字典，每一列是一个基向量
    % y：待稀疏表示的信号
    % K：稀疏度，即选择的基向量个数
    %
    % 输出参数：
    % x_hat：稀疏表示结果
    [~,time] = size(y); %多个信号
    [M, N] = size(D);
    residual_start = zeros(M,time);
    residual_end = zeros(M,time);
    x_hat = zeros(N, time); % 初始化稀疏表示
    for l = 1:time %输入数据有多少列
        x_hat_cnt = zeros(N,1); % 初始化稀疏表示
        y_cnt = y(:,l);
        residual = y_cnt; % 初始化残差
        residual_start(:,l) = residual;
        omega = []; % 选择的基向量索引
        for k = 1:K
            proj = abs(D' * residual); % 计算投影
            [~, idx] = max(proj); % 选择投影最大的基向量
            omega = [omega, idx]; % 更新基向量索引
%             D(:,idx) = zeros(M,1);%清零A的这一列，其实此行可以不要，因为它与残差正交  
            % 利用最小二乘法求解稀疏表示
            x_hat_cnt(omega) = pinv(D(:, omega)) * y_cnt;
    
            % 更新残差
            residual = y_cnt - D(:, omega) * x_hat_cnt(omega);
            
        end
        residual_end(:,l) = residual;
        x_hat(:,l) = x_hat_cnt;
    end
end
% function [ theta ] = omp( y,A,t )  
% % 实现压缩感知OMP算法
% %   y = Phi * x  
% %   x = Psi * theta  
% %   y = Phi*Psi * theta  
% %   t 稀疏度
% %   令 A = Phi*Psi, 则y=A*theta  
% %   现在已知y和A，求theta  
%     [y_rows,y_columns] = size(y);  
%     if y_rows<y_columns  
%         y = y.';%y should be a column vector  
%     end  
%     [M,N] = size(A);%传感矩阵A为M*N矩阵  
%     theta = zeros(N,1);%用来存储恢复的theta(列向量)  
%     At = zeros(M,t);%用来迭代过程中存储A被选择的列  
%     Pos_theta = zeros(1,t);%用来迭代过程中存储A被选择的列序号  
%     r_n = y;%初始化残差(residual)为y
% %     while norm(r_n,"inf")>error_gate
%     for ii=1:t%迭代t次，t为输入参数 
% %         ii = 1;
%         product = A'*r_n;%传感矩阵A各列与残差的内积  
%         [val,pos] = max(abs(product));%找到最大内积绝对值，即与残差最相关的列  
%         At(:,ii) = A(:,pos);%存储这一列  
%         Pos_theta(ii) = pos;%存储这一列的序号  
%         A(:,pos) = zeros(M,1);%清零A的这一列，其实此行可以不要，因为它与残差正交  
%         %y=At(:,1:ii)*theta，以下求theta的最小二乘解(Least Square)  
%         theta_ls = (At(:,1:ii)'*At(:,1:ii))^(-1)*At(:,1:ii)'*y;%最小二乘解  
%         %At(:,1:ii)*theta_ls是y在At(:,1:ii)列空间上的正交投影  
%         r_n = y - At(:,1:ii)*theta_ls;%更新残差
%         norm(r_n,"inf")
%     end  
%     theta(Pos_theta)=theta_ls;%恢复出的theta  
% end  

