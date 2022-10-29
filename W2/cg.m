function [X,error,i] = cg(A,B,X,epsilon,maxit)
% 共轭梯度迭代法
% 输入    - A是一个n*n矩阵
%         - B是一个n*1列向量
%         - X初始值
%         - epsilon允许误差
%         - maxit最大迭代次数
% 输出    - X是方程AX=B的解,error最终误差，i为实际迭代次数
r = B-A*X;
p = r;
for i = 1:maxit
    alpha = (r'*r)/(p'*A*p);
    X = X+alpha*p;
    r1 = r-alpha*A*p;
    beta = (r1'*r1)/(r'*r);
    r = r1;
    p = r+beta*p;
    
    error = norm(A*X-B);
    if error<epsilon
        break;
    end
end
end