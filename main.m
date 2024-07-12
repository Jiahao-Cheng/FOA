%% ��ջ�������
clc;
clear;

%% �������
L = 50;                  % ����߳�
n = 35;                  % �ڵ����
rs = 5;                  % ��֪�뾶
data = 0.4;              % ��ɢ����

%% FOA����
maxgen = 600;          % ��������
sizepop = 20;          % ��Ⱥ��ģ
s = 0.3;               % ����
E0 = 20;             %�ڵ��ʼ����
delta = 200*0.001;   %��λ�����ƶ��ܺ�
beta = 4*0.001;     %�ڵ����smellֵ���ĵ�����
Ec = 20*0.001;       %ÿ�ε����ڵ㸲����������
threshold = 0.0285;      %��ֵ
sigma=2;
E(1,:)=E0+normrnd(0,sigma,[1,n]);
% E = ones(maxgen+1,n) * E0;
% Eo = ones(maxgen+1,n) * E0;

%% �Ľ�FOA
%�����ʼ����ӬȺ��λ��
X_axis = L*rand(1, n);
X0_axis = X_axis;
Y_axis = L*rand(1, n);
Y0_axis = Y_axis;
mov_dist = zeros(maxgen+1,n); %��ʼ���ڵ��ƶ�����
mov_dist0 = zeros(maxgen+1,n); %�Ľ�ǰ��ʼ���ڵ��ƶ�����

% ������ٶ�������Сֵ
for i = 1:sizepop
    % ���λ��
    X(i, :) = X_axis + 2*s*rand(1, n)-s;
    Y(i, :) = Y_axis + 2*s*rand(1, n)-s;
    % ζ��Ũ�Ⱥ���(������)
    Smell(i) = computeSmell(X(i, :), Y(i, :), L, rs, data);
end
% �ҳ��˹�ӬȺ����ζ��Ũ����ߵĹ�Ӭ(�󼫴�ֵ)
[Smell_new, Index_new]=max(Smell);
% �����ζŨ�ȡ���Ӭλ�á���Ӧ������λ��
X_axis = X(Index_new, :);
Y_axis = Y(Index_new, :);
Smellbest = Smell_new;
E0=E0-beta-Ec;
% E = ones(maxgen+1,n) * E0;
T_life(1) = min(E(1,:)) / Ec;
P(1)=Smellbest * T_life(1);

% ��ʼ�����ʾ
gbest = [X_axis; Y_axis]';
disp('��ʼλ�ã�' );
disp([num2str(gbest)]);
disp(['��ʼ�����ʣ�', num2str(Smellbest)]);
% ��ʼ����ͼ
figure;
for i = 1:n
    axis([0 L 0 L]);            % �������귶Χ
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % �Ƕ�[0, 2*pi]
    hold on;
    p2 = fill(x(i)+rs*cos(sita), y(i)+rs*sin(sita), 'y');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSN�ڵ�', '��������'});
title 'FOA-WSN��ʼ���';

%% �Ľ�FOA��Ӭ����Ѱ��
for gen = 2:(maxgen+1)
    %�жϽڵ㲢ɸ�������ڵ�
    %��ʼ���ڵ�����
%     E(gen,:)=E(gen-1,:);
    % ����λ�ú��ٶȸ���
    for i = 1:sizepop
        X(i, :) = X_axis + 2*s*rand(1, n)-s;
        Y(i, :) = Y_axis + 2*s*rand(1, n)-s;
        % �߽紦��
        X(i, :) = max(X(i, :), 0);
        X(i, :) = min(X(i, :), L);
        Y(i, :) = max(Y(i, :), 0);
        Y(i, :) = min(Y(i, :), L);
        % ���㸲����
        Smell(i) = computeSmell(X(i, :), Y(i, :), L, rs, data);
    end
    % ������ζŨ��ֵѰ�Ҽ�ֵ
    [Smell_new, Index_new]=max(Smell);
    % �������ֵλ��
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
    % ����ÿ�ε�����ÿ���ڵ�ʣ����������������ʱ����ۺ�����
            E_mov(gen,i) = delta * mov_dist(gen,i); %�ڵ��ƶ����ĵ�����
            E_compute =gen * beta; %�ڵ����smell���ĵ�������
            E_cover=gen * Ec; %�ڵ㸲�����ĵ�������
            E(gen,i) = E(1,i) - E_compute- E_cover -sum(E_mov(:,i));%�ڵ���j�ֵ����������
        end
    end
    % ÿ������Smellֵ��¼��yy�����У�����¼���ŵ�������
    yy(gen-1) = Smellbest;
    Xbest(gen-1, :) = X_axis;
    Ybest(gen-1, :) = Y_axis;
    % ��ʾ������Ϣ
    display(['FOA: At iteration ', num2str(gen-1), ' the best fitness is ', num2str(yy(gen-1))]);
    % ����˴ε��������������ں��ۺ�����
    for i=1:n
        if E(gen,i) < 0
            E(gen,i)=0;
        end    
    end    
    T_life(gen) = min(E(gen,:)) / Ec;%ÿ�ε��������������ʱ��
    P(gen)=Smellbest * T_life(gen);%�ۺ�����
    if P(gen) <= 0
        break
    end    
end

% �����ʾ
gbest = [ Xbest(end, :); Ybest(end, :)]';
disp('����λ�ã�');
disp([num2str(gbest)]);
disp(['���Ÿ����ʣ�', num2str(yy(end))]);

%% FOA
%������ٶ�������Сֵ
for i = 1:sizepop
    % ���λ��
    X0(i, :) = X0_axis + 2*s*rand(1, n)-s;
    Y0(i, :) = Y0_axis + 2*s*rand(1, n)-s;
    % ζ��Ũ�Ⱥ���(������)
    Smell0(i) = computeSmell(X0(i, :), Y0(i, :), L, rs, data);
end
% �ҳ��˹�ӬȺ����ζ��Ũ����ߵĹ�Ӭ(�󼫴�ֵ)
[Smell_new, Index_new]=max(Smell0);
% �����ζŨ�ȡ���Ӭλ�á���Ӧ������λ��
X0_axis = X0(Index_new, :);
Y0_axis = Y0(Index_new, :);
Smellbest0 = Smell_new;
% Eo = ones(maxgen+1,n) * E0;
Eo(1,:)=E(1,:);
T_life0(1) = min(Eo(1,:)) / Ec;
P0(1)=Smellbest0 * T_life0(1);

% ��ʼ�����ʾ
gbest0 = [X0_axis; Y0_axis]';
disp('��ʼλ�ã�' );
disp([num2str(gbest0)]);
disp(['��ʼ�����ʣ�', num2str(Smellbest0)]);


%% FOA��Ӭ����Ѱ��
for gen = 2:(maxgen+1)
    %�жϽڵ㲢ɸ�������ڵ�
%     Eo(gen,:)=Eo(gen-1,:);
    % ����λ�ú��ٶȸ���
    for i = 1:sizepop
        X0(i, :) = X0_axis + 2*s*rand(1, n)-s;
        Y0(i, :) = Y0_axis + 2*s*rand(1, n)-s;
        % �߽紦��
        X0(i, :) = max(X0(i, :), 0);
        X0(i, :) = min(X0(i, :), L);
        Y0(i, :) = max(Y0(i, :), 0);
        Y0(i, :) = min(Y0(i, :), L);
        % ���㸲����
        Smell0(i) = computeSmell(X0(i, :), Y0(i, :), L, rs, data);
    end
    % ������ζŨ��ֵѰ�Ҽ�ֵ
    [Smell_new, Index_new]=max(Smell0);
    % �������ֵλ��
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
    % ����ÿ�ε�����ÿ���ڵ�ʣ����������������ʱ����ۺ�����
            Eo_mov(gen,i) = delta * mov_dist0(gen,i); %�ڵ��ƶ����ĵ�����
            Eo_compute =gen * beta; %�ڵ����smell���ĵ�������
            Eo_cover=gen * Ec; %�ڵ㸲�����ĵ�������
            Eo(gen,i) = Eo(1,i) - Eo_compute- Eo_cover -sum(Eo_mov(:,i));%�ڵ���j�ֵ����������
        end
    end
    % ÿ������Smellֵ��¼��yy�����У�����¼���ŵ�������
    yy0(gen-1) = Smellbest0;
    X0best(gen-1, :) = X0_axis;
    Y0best(gen-1, :) = Y0_axis;
    % ��ʾ������Ϣ
    display(['FOA0: At iteration ', num2str(gen-1), ' the best fitness is ', num2str(yy0(gen-1))]);
    % ����˴ε������нڵ�����
    for i=1:n
        if Eo(gen,i) < 0
            Eo(gen,i)=0;
        end    
    end    
    T_life0(gen) = min(Eo(gen,:)) / Ec;%ÿ�ε��������������ʱ��
    P0(gen)=Smellbest0 * T_life0(gen);%�ۺ�����
    if P0(gen) <= 0
        break
    end    
end

% �����ʾ
gbest0 = [ X0best(end, :); Y0best(end, :)]';
disp('����λ�ã�');
disp([num2str(gbest0)]);
disp(['���Ÿ����ʣ�', num2str(yy0(end))]);

%% ��ͼ
figure;
plot(yy, 'r', 'lineWidth', 2);  %  ��������ͼ
hold on;
plot(yy0, 'b', 'lineWidth', 2);
title('������ѵ������', 'fontsize', 12);
xlabel('��������', 'fontsize', 12);
ylabel('������', 'fontsize', 12);
legend('�Ľ�FOA','FOA');

figure;
plot(P, 'r', 'lineWidth', 2);
hold on;
plot(P0, 'b', 'lineWidth', 2);
title('�ۺ�����Pѵ������', 'fontsize', 12);
xlabel('��������', 'fontsize', 12);
ylabel('�ۺ�����P', 'fontsize', 12);
legend('�Ľ�FOA','FOA');

figure;
plot(T_life, 'r', 'lineWidth', 2);
hold on;
plot(T_life0, 'b', 'lineWidth', 2);
title('��������ѵ������', 'fontsize', 12);
xlabel('��������', 'fontsize', 12);
ylabel('��������T_life', 'fontsize', 12);
legend('�Ľ�FOA','FOA');

figure;
for i = 1:n
    axis([0 L 0 L]);            % �������귶Χ
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % �Ƕ�[0, 2*pi]
    hold on;
    p2 = fill(x(i)+rs*cos(sita), y(i)+rs*sin(sita), 'g');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSN�ڵ�', '��������'});
title '�Ľ�FOA�㷨���ս��';
figure;
for i = 1:n
    axis([0 L 0 L]);            % �������귶Χ
    x0 = gbest0(:, 1);
    y0 = gbest0(:, 2);
    sita = 0:pi/100:2*pi;   % �Ƕ�[0, 2*pi]
    hold on;
    p2 = fill(x0(i)+rs*cos(sita), y0(i)+rs*sin(sita), 'g');
end
p1 = plot(gbest0(:, 1), gbest0(:, 2), 'r*');
legend([p1, p2], {'WSN�ڵ�', '��������'});
title 'FOA�㷨���ս��';