function [P,err,iter] = newdim(F,JF,P,delta,epsilon,maxit)
% 牛顿-拉夫森法求解非线性方程组F(X)=0,给定初始近似值P0，收敛到解P
% 输入    - F非线性系统
%         - JF为非线性系统的雅可比矩阵
%         - P为解的初始近似值
%         - delta和epsilon分别为P和F(P)的允许误差
%         - maxit为最大迭代次数
% 输出    - P为解的近似值
%         - err为P的估计误差，iter为实际迭代次数
Y = F(P);
for iter = 1:maxit
    J = JF(P);
    Q = P-(J\Y')';
    err = norm(Q-P);
    relerr = err/(norm(Q)+eps);
    P = Q;
    Y = F(P);
    if (err<delta) || (relerr<delta) || (norm(Y)<epsilon)
        break;
    end
end
end