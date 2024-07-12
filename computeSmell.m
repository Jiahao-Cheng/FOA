function z = computeSmell(x, y, L, R, data)
%% 适应度函数：WSN的覆盖率
% input：
% x        圆心横坐标
% y        圆心纵坐标
% L        区域边长
% R        通信半径
% data     离散粒度
% output:
% z        覆盖率
N = length(x);                      % 节点总个数
[m, n] = meshgrid(0:data:L);        % 离散化区域内的点
[row, col] = size(m);
for i = 1:N
    D = sqrt((m-x(i)).^2+(n-y(i)).^2);     % 计算坐标点到圆心的距离
    [m0, n0] = find(D <= R);               % 检测出圆覆盖点的坐标
    Ind = (m0-1).*col+n0;                  % 坐标与索引转化
    M(Ind) = 1;                            % 改变覆盖状态
end
scale = sum(M(1:end))/(row*col);           % 计算覆盖比例
z = scale;

