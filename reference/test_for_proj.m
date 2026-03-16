clear,clc
dx = 0.5; % x方向步长
dy = 0.5; % y方向步长
sx = 5;  % x方向抽样间隔
sy = 5;  % y方向抽样间隔
%% 定义
% 约束定义：
theta = 10;     % 临界半径
epsilon = 10;    % 极限膨胀层
f = @(x) (x(1).^2 + x(2).^2 - theta^2)/(2*theta*epsilon); % 约束函数
% 势场函数：
x_m = 25;
y_m = 10;
phi = @(x) (x(1)-x_m).^2 + 2*(x(2)-y_m).^2;
syms x [2,1];
grad_phi = jacobian(phi(x),x);
grad_handle = matlabFunction(grad_phi,'Vars', {x});
grad_phi = @(x)grad_handle(x)';
% 计算区域
X = -25:dx:25;
Y = -25:dy:25;
Z = zeros(length(Y),length(X));
grad_Z = zeros(length(Y),length(X),2);
grad_P = zeros(length(Y),length(X),2);
constr = zeros(length(Y),length(X));
%% 计算
for i = 1:length(Y)
    tic
    for j = 1:length(X)
        x = X(j);
        y = Y(i);
        % 计算势函数：
        Z(i,j) = phi([x;y]);
        % 计算梯度：
        grad_Z(i,j,:) = grad_phi([x;y]);
        % 计算投影结果：
        grad_P(i,j,:) = Proj([x;y],grad_phi([x;y]),f);
        % 计算约束范围：
        if f([x;y])<0   % 约束内
            constr(i,j) = -1;
        elseif f([x;y])<1   % 临界约束
            constr(i,j) = 0;
        else    % 约束外
            constr(i,j) = 1;
        end
    end
    disp(toc)
end


%% 3D可视化
% 默认解释器 - latex
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
figure(1)
[x_grid, y_grid] = meshgrid(X,Y);

% 势函数：
surf(x_grid,y_grid,Z,'FaceAlpha',0.8,'EdgeColor','none')

hold on 
plot3(x_m,y_m,0,'rx','LineWidth',2)

% 梯度：
indx = 1:sx:length(X);
indy = 1:sy:length(Y);
x_quiv = x_grid(indx,indx);
y_quiv = y_grid(indy,indy);
z_quiv = Z(indx,indy);
quiver3(x_quiv, y_quiv, z_quiv, grad_Z(indx,indy,1), grad_Z(indx,indy,2), zeros(length(indx),length(indy)),...
        0.02, 'LineWidth', 1.5, 'Color', 'r', 'MaxHeadSize', 0.5)
% 装饰：
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
%% 2D等高线
figure(2)
% 势函数:
contour(x_grid,y_grid,Z,20);
hold on
plot(x_m,y_m,'rx','LineWidth',2,'MarkerSize',10)
% 约束：
[C, h] = contour(x_grid,y_grid,constr,[0, 1],'LineStyle','--','ShowText','on','LabelFormat','$f(x)=%.0f$');
clabel(C, h, [0, 1], 'Interpreter', 'latex', 'FontSize', 16, 'LabelSpacing', inf);
% 梯度：
quiver(x_quiv,y_quiv,grad_Z(indx,indy,1),grad_Z(indx,indy,2),...
       1, 'LineWidth', 1.5, 'Color', 'r', 'MaxHeadSize', 0.5)
% 投影梯度：
quiver(x_quiv,y_quiv,grad_P(indx,indy,1),grad_P(indx,indy,2),...
       4, 'LineWidth', 1.5, 'Color', 'b', 'MaxHeadSize', 0.5)
% 装饰：
xlabel('$x$')
ylabel('$y$')
axis equal
