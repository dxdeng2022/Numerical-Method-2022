function [p0,err,i,y,p_mat] = newton(f,df,p0,delta,epsilon,maxi)
% 牛顿迭代法
% 输入    - f,df分别为目标函数及其导数
%         - delta,epsilon分别为p0、y的允许误差
%         - maxi为最大迭代次数
% 输出    - p0为迭代结果
%         - err 为p0的误差值
%         - i迭代次数，y为迭代结果对应的y值
p_mat = zeros(1,maxi);
p_mat(1) = p0;
for i = 1:maxi
    p1 = p0 - f(p0)/df(p0);
    err = abs(p1 - p0);
    relaerr = 2*err/(abs(p1)+delta);
    p0 = p1;
    y = f(p0);
    p_mat(i+1) = p0;
    if err <= delta || relaerr <= delta || abs(y) <= epsilon
        break;
    end
end
end