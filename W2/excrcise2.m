% 计算超松弛因子  -- 2022-10-17
n = 100;
e = ones(n,1);
A = spdiags([e -2*e 10*e -2*e e],-2:2,n,n);
D = diag(diag(A));     % 对角元素   
L = tril(A,-1);        % 下三角元素
U = triu(A,1);         % 上三角元素
w = 1:0.01:2;          % 超松弛因子omiga
P_w = zeros(size(w));  % 谱半径列表
for i = 1:length(w)
    tmp = D+w(i)*L;
    B_w = tmp\((1-w(i))*D-w(i)*U);
    p_w = max(abs(eig(B_w)));      % 计算谱半径
    P_w(i) = p_w;
end
figure; box on; grid on; hold on;
plot(w,P_w,'k','LineWidth',1.5);
xlabel("\omega"); ylabel("\rho (\omega)");
set(gca,'fontname','Times New Roman','fontsize',14);
w_opt = w(P_w == min(P_w));        % 最优超松弛因子
scatter(w_opt,min(P_w),'LineWidth',1.5);

%% % excrcise 2-1(a) -- 2022-10-20
m = 1:100;
time_data = zeros(5,length(m));           % 迭代时间数据
iter_data = zeros(5,length(m));           % 迭代次数数据
for i = 1:length(m)
n = 100*m(i);
e = ones(n,1);
A = spdiags([e -2*e 10*e -2*e e],-2:2,n,n);
b = repmat([-3;3], [n/2 1]);
% D = diag(diag(A));              % 对角元素
% L = tril(A,-1);                 % 下三角元素
% pre_mat = (D+0.5*L)/sqrtm(D)/sqrt(0.5*(2-0.5));    % pcg预处理矩阵
relerr = 5e-4/norm(b);          % 相对残差

if m(i) <= 20
    tic;
    X = uptrbk(A,b);                               % 高斯消去法,m=25需100s左右
    time_data(1,i) = toc;
end
X0 = ones(n,1);                                    % 初始值为1
tic
[~,~,iter1] = sor(A,b,1.03,X0,5e-4,1000);          % SOR方法
time_data(2,i) = toc; tic;
[~,~,iter2] = cg(A,b,X0,5e-4,1000);                % 共轭梯度法
time_data(3,i) = toc; tic;
[~,~,~,iter3] = pcg(A,b,relerr,1000,[],[],X0);        % pcg
time_data(4,i) = toc; tic
[~,~,~,iter4] = gmres(A,b,[],relerr,100,[],[],X0);    % gmres
time_data(5,i) = toc;
iter_data(:,i) = [n; iter1; iter2; iter3; iter4(1)*iter4(2)];
end

% 迭代时间绘图
% 最大图
figure; box on;
plot(m(1:20),time_data(1,1:20),'k--','linewidth',1.5);    % 高斯消去耗时曲线
set(gca,'fontname','Times New Roman','fontsize',14);
legend("Gauss",'Location','southeast'); legend('boxoff');
xlabel("m"); ylabel("Times (s)");

% 小图
axes('Position',[0.1764,0.2476,0.6336,0.6568]); box on; hold on;
plot(m,time_data(2,:),'k','linewidth',1.5);               % SOR耗时曲线
plot(m,time_data(3,:),'r','linewidth',1.5);               % cg耗时曲线
plot(m,time_data(4,:),'g','linewidth',1.5);               % pcg耗时曲线
plot(m,time_data(5,:),'b','linewidth',1.2);               % gmres耗时曲线
legend(["SOR","cg","pcg","gmres"]);
legend("boxoff");
set(gca,'fontname','Times New Roman','fontsize',14);

% 小图局部放大
axes('Position',[0.25,0.5771,0.3629,0.2981]); box on; hold on;
plot(m,time_data(4,:),'g','linewidth',1.5);
plot(m,time_data(5,:),'b','linewidth',1.2);
set(gca,'fontname','Times New Roman','fontsize',12);

% 迭代次数绘图
figure; box on; hold on;
plot(m,iter_data(1,:),'k--','linewidth',1.5);         % 高斯消去迭代次数曲线
legend("Gauss",'location','southeast');
legend("boxoff");
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel("m"); ylabel("IterationNumber");

% 小图
axes('Position',[0.20,0.31,0.52,0.56]); box on; hold on;
plot(m,iter_data(2,:),'k','linewidth',1.5);           % SOR迭代次数曲线
plot(m,iter_data(3,:),'-ro','linewidth',2);           % cg迭代次数曲线
plot(m,iter_data(4,:),'g','linewidth',1.5);           % pcg迭代次数曲线
plot(m,iter_data(5,:),'--b','linewidth',0.8);             % gmres迭代次数曲线
legend(["SOR","cg","pcg","gmres"]);
legend("boxoff");
set(gca,'fontname','Times New Roman','fontsize',14);


%% exercise 2-1(b) -- 2022-10-28
clear; close all;
m = [1 5 10 50 100];
err_data = zeros(4*length(m),100);     % 残差序列,4种方法，每4行存m的一种情况
iter_data = zeros(size(err_data,1),1); % 迭代次数数据
for i = 1:length(m)
    n = 100*m(i);
    e = ones(n,1);
    A = spdiags([e -2*e 10*e -2*e e],-2:2,n,n);
    b = repmat([-3;3], [n/2 1]);
    D = diag(diag(A));                  % 对角元素
    L = -tril(A,-1);                    % 负下三角元素
    D_sqrt = sqrtm(D);                  % 对角阵开根号
    w = 0.5;
    pre_mat = (D-w*L)/D_sqrt/sqrt(w*(2-w));   % 预处理矩阵
    relerr = 5e-4/norm(b);

    X0 = ones(n,1);
    [~,~,iter1,err1] = sor(A,b,1.03,X0,5e-4,100);                     % sor
    [~,~,iter2,err2] = cg(A,b,X0,5e-4,100);                           % cg
    [~,~,~,iter3,err3] = pcg(A,b,relerr,100,pre_mat,pre_mat',X0);     % pcg
    [~,~,~,iter4,err4] = gmres(A,b,[],relerr,100,pre_mat,eye(n),X0);  % gmres
    % 存储残差数据和迭代次数
    err_data(4*i-3,1:iter1+1) = err1;
    err_data(4*i-2,1:iter2+1) = err2;
    err_data(4*i-1,1:iter3+1) = err3;
    err_data(4*i,1:iter4(1)*iter4(2)+1) = err4;
    iter_data(4*i-3:4*i) = [iter1; iter2; iter3; iter4(1)*iter4(2)];
end

% 绘图程序部分
for j = 1:length(m)
    figure; box on; hold on;
    tmp = 1:iter_data(4*j-3)+1;
    plot(tmp-1,err_data(4*j-3,tmp),'-ko','LineWidth',1.5);    % sor残差
    tmp = 1:iter_data(4*j-2)+1;
    plot(tmp-1,err_data(4*j-2,tmp),'-ro','LineWidth',1.5);    % cg残差
    tmp = 1:iter_data(4*j-1)+1;
    plot(tmp-1,err_data(4*j-1,tmp),'-go','LineWidth',1);      % pcg残差
    tmp = 1:iter_data(4*j)+1;
    plot(tmp-1,err_data(4*j,tmp),'-bo','LineWidth',1.5);      % gmres残差
    set(gca,'fontname','Times New Roman','FontSize',18);
    legend(["sor","cg","pcg","gmres"]);
    legend("boxoff");
    xlabel("k"); ylabel("||AX-b||");
    title(sprintf("m = %d",m(j)));
end


%% excrcise 2-1(c) -- 2022-10-23
clear; close all;
m = [1 10 100];
err_data = zeros(2*length(m),100);        % 残差序列
iter_data = zeros(size(err_data,1),1);    % 迭代次数
for i = 1:length(m)
    n = 100*m(i);
    e = ones(n,1);
    A = spdiags([e -2*e 10*e -2*e e],-2:2,n,n);
    b = repmat([-3;3], [n/2 1]);
    D = diag(diag(A));                  % 对角元素
    L = -tril(A,-1);                    % 负下三角元素
    D_sqrt = sqrtm(D);                  % 对角阵开根号
    w = 0.5;
    C = (D-w*L)/D_sqrt/sqrt(w*(2-w));   % 预处理矩阵
    A_pre = C\A/(C');                   % 即inv(C)*A*inv(C');
    A_cond = cond(A);
    A_pre_cond = cond(A_pre);
    fprintf("m = %d: cond(A)=%.4f, cond(A_pre)=%.4f\n",m(i),A_cond,A_pre_cond);
    
    X0 = ones(n,1);
    [~,~,iter1,err1] = cg(A,b,X0,5e-4,100);
    [~,~,~,iter2,err2] = pcg(A,b,5e-4/norm(b),100,C,C',X0);
    err_data(2*i-1,1:iter1+1) = err1;
    err_data(2*i,1:iter2+1) = err2;
    iter_data(2*i-1:2*i) = [iter1; iter2];
end

% 绘图程序
figure; box on; hold on;
for j = 2:3
    tmp = 1:iter_data(2*j-1)+1;
    plot(tmp-1,err_data(2*j-1,tmp),'-ko','LineWidth',1.5);
    tmp = 1:iter_data(2*j)+1;
    plot(tmp-1,err_data(2*j,tmp),'-ro','LineWidth',1.2);
end
legend(["cg","pcg"]);legend("boxoff");
xlabel("k"); ylabel("||Ax-b||");
set(gca,'fontname','Times New Roman','fontsize',14);

% 小图
axes("position",[0.3414,0.2562,0.5272,0.5553]); box on; hold on;
for j = 1
    tmp = 1:iter_data(2*j-1);
    plot(tmp-1,err_data(2*j-1,tmp),'-ko','LineWidth',1.5);
    tmp = 1:iter_data(2*j);
    plot(tmp-1,err_data(2*j,tmp),'-ro','LineWidth',1.5);
end
set(gca,'fontname','Times New Roman','fontsize',14);



%% % excrcise 2-2 -- 2022-10-18
x1 = -1.5:0.01:1.5;
y1 = 7*x1.^3-10*x1-1;
y2 = x1;
x2 = -8*y2.^3+11*y2+1;
figure; box on; hold on;
plot(x1,y1,'k','LineWidth',1.5);
plot(x2,y2,'r','LineWidth',1.5);
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel("x"); ylabel("y");
set(gca,'XMinorTick','on','YMinorTick','on');

% 通过绘图得到的交点
intersect_x = [-1.06, -0.231, 1.291, -1.153, -0.091, 1.243, -1.198, 0.0125, 1.1861];
intersect_y = [1.257, 1.225, 1.159, -0.202, -0.1, 0.0221, -1.0558, -1.1248, -1.181];
scatter(intersect_x,intersect_y,'*k');

% 牛顿-拉夫森迭代
F = @func_W2_2;
JF = @func_W2_2_df;
intersect_points = zeros(2,9);
for i = 1:9
    % 添加0.1以内的随机扰动
    P0 = [intersect_x(i)+unifrnd(-0.1,0.1),intersect_y(i)+unifrnd(-0.1,0.1)];    
    [P,err,iter] = newdim(F,JF,P0,5e-9,5e-9,1000);
    intersect_points(:,i) = P';
end


% 非线性方程组  --  2022-10-18
function Z = func_W2_2(X)
x = X(1); y = X(2);
Z = zeros(1,2);
Z(1) = 7*x.^3-10*x-y-1;
Z(2) = 8*y.^3-11*y+x-1;
end

% 非线性方程组的雅克比矩阵  -- 2022-10-18
function JZ = func_W2_2_df(X)
x = X(1); y = X(2);
JZ = [21*x.^2-10 -1; 1 24*y.^2-11];
end