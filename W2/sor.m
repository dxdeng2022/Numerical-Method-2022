function [X,error,i] = sor(A,B,omiga,X,epsilon,maxit)
% 超松弛迭代算法
% 输入    - A是一个n*n矩阵
%         - B是一个n*1列向量
%         - omiga超松弛因子
%         - X初始值
%         - epsilon允许误差
%         - maxit最大迭代次数
% 输出    - X是方程AX=B的解，error最终误差，i为实际迭代次数
D = diag(diag(A));     % 对角元素   
L = tril(A,-1);        % 下三角元素
U = triu(A,1);         % 上三角元素
M = (D+omiga*L)\((1-omiga)*D-omiga*U);
beta = (D+omiga*L)\(omiga*B);
for i = 1:maxit
   X = M*X+beta;
   error = norm(A*X-B);
   if error < epsilon
       break;
   end
end
end