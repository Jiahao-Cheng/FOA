%% 清空环境变量
clc;
clear;

%% 网络参数
L = 50;                  % 区域边长
n = 35;                  % 节点个数
rs = 5;                  % 感知半径
data = 0.4;              % 离散粒度

%% FOA参数
maxgen = 600;          % 迭代次数
sizepop = 20;          % 种群规模
s = 0.3;               % 步长
E0 = 20;             %节点初始能量
delta = 200*0.001;   %单位距离移动能耗
beta = 4*0.001;     %节点计算smell值消耗的能量
Ec = 20*0.001;       %每次迭代节点覆盖所需能量
threshold = 0.0285;      %阈值
sigma=2;
E(1,:)=E0+normrnd(0,sigma,[1,n]);
% E = ones(maxgen+1,n) * E0;
% Eo = ones(maxgen+1,n) * E0;

%% 改进FOA
%随机初始化果蝇群体位置
X_axis = L*rand(1, n);
X0_axis = X_axis;
Y_axis = L*rand(1, n);
Y0_axis = Y_axis;
mov_dist = zeros(maxgen+1,n); %初始化节点移动距离
mov_dist0 = zeros(maxgen+1,n); %改进前初始化节点移动距离

% 个体和速度最大和最小值
for i = 1:sizepop
    % 随机位置
    X(i, :) = X_axis + 2*s*rand(1, n)-s;
    Y(i, :) = Y_axis + 2*s*rand(1, n)-s;
    % 味道浓度函数(覆盖率)
    Smell(i) = computeSmell(X(i, :), Y(i, :), L, rs, data);
end
% 找出此果蝇群体中味道浓度最高的果蝇(求极大值)
[Smell_new, Index_new]=max(Smell);
% 最佳气味浓度、果蝇位置、适应度最优位置
X_axis = X(Index_new, :);
Y_axis = Y(Index_new, :);
Smellbest = Smell_new;
E0=E0-beta-Ec;
% E = ones(maxgen+1,n) * E0;
T_life(1) = min(E(1,:)) / Ec;
P(1)=Smellbest * T_life(1);

% 初始结果显示
gbest = [X_axis; Y_axis]';
disp('初始位置：' );
disp([num2str(gbest)]);
disp(['初始覆盖率：', num2str(Smellbest)]);
% 初始覆盖图
figure;
for i = 1:n
    axis([0 L 0 L]);            % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i)+rs*cos(sita), y(i)+rs*sin(sita), 'y');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSN节点', '覆盖区域'});
title 'FOA-WSN初始结果';

%% 改进FOA果蝇迭代寻优
for gen = 2:(maxgen+1)
    %判断节点并筛除死亡节点
    %初始化节点能量
%     E(gen,:)=E(gen-1,:);
    % 粒子位置和速度更新
    for i = 1:sizepop
        X(i, :) = X_axis + 2*s*rand(1, n)-s;
        Y(i, :) = Y_axis + 2*s*rand(1, n)-s;
        % 边界处理
        X(i, :) = max(X(i, :), 0);
        X(i, :) = min(X(i, :), L);
        Y(i, :) = max(Y(i, :), 0);
        Y(i, :) = min(Y(i, :), L);
        % 计算覆盖率
        Smell(i) = computeSmell(X(i, :), Y(i, :), L, rs, data);
    end
    % 根据气味浓度值寻找极值
    [Smell_new, Index_new]=max(Smell);
    % 保留最佳值位置
    if Smell_new > Smellbest
        for i = 1:n
            if E(gen-1,i) > 0
              if ((E(gen-1,i)/sum(E(gen-1,:))) >= threshold) 
                  mov_dist(gen,i) = sqrt((X_axis(i)-X(Index_new,i)).^2+(Y_axis(i)-Y(Index_new,i)).^2);
                  X_axis(i) = X(Index_new, i);
                  Y_axis(i) = Y(Index_new, i);
                  Smellbest = Smell_new;         
              end
            else
                E(gen,i) = 0;
            end
        end
    end
    for i=1:n
        if E(gen-1,i) > 0
    % 计算每次迭代后每个节点剩余能量，网络生存时间和综合性能
            E_mov(gen,i) = delta * mov_dist(gen,i); %节点移动消耗的能量
            E_compute =gen * beta; %节点计算smell消耗的总能量
            E_cover=gen * Ec; %节点覆盖消耗的总能量
            E(gen,i) = E(1,i) - E_compute- E_cover -sum(E_mov(:,i));%节点在j轮迭代后的能量
        end
    end
    % 每代最优Smell值记录到yy数组中，并记录最优迭代坐标
    yy(gen-1) = Smellbest;
    Xbest(gen-1, :) = X_axis;
    Ybest(gen-1, :) = Y_axis;
    % 显示迭代信息
    display(['FOA: At iteration ', num2str(gen-1), ' the best fitness is ', num2str(yy(gen-1))]);
    % 计算此次迭代网络生命周期和综合性能
    for i=1:n
        if E(gen,i) < 0
            E(gen,i)=0;
        end    
    end    
    T_life(gen) = min(E(gen,:)) / Ec;%每次迭代后的网络生命时间
    P(gen)=Smellbest * T_life(gen);%综合性能
    if P(gen) <= 0
        break
    end    
end

% 结果显示
gbest = [ Xbest(end, :); Ybest(end, :)]';
disp('最优位置：');
disp([num2str(gbest)]);
disp(['最优覆盖率：', num2str(yy(end))]);

%% FOA
%个体和速度最大和最小值
for i = 1:sizepop
    % 随机位置
    X0(i, :) = X0_axis + 2*s*rand(1, n)-s;
    Y0(i, :) = Y0_axis + 2*s*rand(1, n)-s;
    % 味道浓度函数(覆盖率)
    Smell0(i) = computeSmell(X0(i, :), Y0(i, :), L, rs, data);
end
% 找出此果蝇群体中味道浓度最高的果蝇(求极大值)
[Smell_new, Index_new]=max(Smell0);
% 最佳气味浓度、果蝇位置、适应度最优位置
X0_axis = X0(Index_new, :);
Y0_axis = Y0(Index_new, :);
Smellbest0 = Smell_new;
% Eo = ones(maxgen+1,n) * E0;
Eo(1,:)=E(1,:);
T_life0(1) = min(Eo(1,:)) / Ec;
P0(1)=Smellbest0 * T_life0(1);

% 初始结果显示
gbest0 = [X0_axis; Y0_axis]';
disp('初始位置：' );
disp([num2str(gbest0)]);
disp(['初始覆盖率：', num2str(Smellbest0)]);


%% FOA果蝇迭代寻优
for gen = 2:(maxgen+1)
    %判断节点并筛除死亡节点
%     Eo(gen,:)=Eo(gen-1,:);
    % 粒子位置和速度更新
    for i = 1:sizepop
        X0(i, :) = X0_axis + 2*s*rand(1, n)-s;
        Y0(i, :) = Y0_axis + 2*s*rand(1, n)-s;
        % 边界处理
        X0(i, :) = max(X0(i, :), 0);
        X0(i, :) = min(X0(i, :), L);
        Y0(i, :) = max(Y0(i, :), 0);
        Y0(i, :) = min(Y0(i, :), L);
        % 计算覆盖率
        Smell0(i) = computeSmell(X0(i, :), Y0(i, :), L, rs, data);
    end
    % 根据气味浓度值寻找极值
    [Smell_new, Index_new]=max(Smell0);
    % 保留最佳值位置
    if Smell_new > Smellbest0
        for i = 1:n
            if Eo(gen-1,i) > 0
                  mov_dist0(gen,i) = sqrt((X0_axis(i)-X0(Index_new,i)).^2+(Y0_axis(i)-Y0(Index_new,i)).^2);
                  X0_axis(i) = X0(Index_new, i);
                  Y0_axis(i) = Y0(Index_new, i);
                  Smellbest0 = Smell_new;   
            else
                Eo(gen,i) = 0;
            end
        end
    end
    for i=1:n
        if Eo(gen-1,i) > 0
    % 计算每次迭代后每个节点剩余能量，网络生存时间和综合性能
            Eo_mov(gen,i) = delta * mov_dist0(gen,i); %节点移动消耗的能量
            Eo_compute =gen * beta; %节点计算smell消耗的总能量
            Eo_cover=gen * Ec; %节点覆盖消耗的总能量
            Eo(gen,i) = Eo(1,i) - Eo_compute- Eo_cover -sum(Eo_mov(:,i));%节点在j轮迭代后的能量
        end
    end
    % 每代最优Smell值记录到yy数组中，并记录最优迭代坐标
    yy0(gen-1) = Smellbest0;
    X0best(gen-1, :) = X0_axis;
    Y0best(gen-1, :) = Y0_axis;
    % 显示迭代信息
    display(['FOA0: At iteration ', num2str(gen-1), ' the best fitness is ', num2str(yy0(gen-1))]);
    % 计算此次迭代所有节点能量
    for i=1:n
        if Eo(gen,i) < 0
            Eo(gen,i)=0;
        end    
    end    
    T_life0(gen) = min(Eo(gen,:)) / Ec;%每次迭代后的网络生命时间
    P0(gen)=Smellbest0 * T_life0(gen);%综合性能
    if P0(gen) <= 0
        break
    end    
end

% 结果显示
gbest0 = [ X0best(end, :); Y0best(end, :)]';
disp('最优位置：');
disp([num2str(gbest0)]);
disp(['最优覆盖率：', num2str(yy0(end))]);

%% 绘图
figure;
plot(yy, 'r', 'lineWidth', 2);  %  画出迭代图
hold on;
plot(yy0, 'b', 'lineWidth', 2);
title('覆盖率训练过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('覆盖率', 'fontsize', 12);
legend('改进FOA','FOA');

figure;
plot(P, 'r', 'lineWidth', 2);
hold on;
plot(P0, 'b', 'lineWidth', 2);
title('综合性能P训练过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('综合性能P', 'fontsize', 12);
legend('改进FOA','FOA');

figure;
plot(T_life, 'r', 'lineWidth', 2);
hold on;
plot(T_life0, 'b', 'lineWidth', 2);
title('网络生命训练过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('网络生命T_life', 'fontsize', 12);
legend('改进FOA','FOA');

figure;
for i = 1:n
    axis([0 L 0 L]);            % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i)+rs*cos(sita), y(i)+rs*sin(sita), 'g');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSN节点', '覆盖区域'});
title '改进FOA算法最终结果';
figure;
for i = 1:n
    axis([0 L 0 L]);            % 限制坐标范围
    x0 = gbest0(:, 1);
    y0 = gbest0(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x0(i)+rs*cos(sita), y0(i)+rs*sin(sita), 'g');
end
p1 = plot(gbest0(:, 1), gbest0(:, 2), 'r*');
legend([p1, p2], {'WSN节点', '覆盖区域'});
title 'FOA算法最终结果';