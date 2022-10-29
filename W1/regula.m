function [c,err,i,yc,c_mat] = regula(f,a,b,delta,epsilon,maxi)
% 试值迭代法
% 输入    - f为目标迭代函数
%         - [a, b]为初始迭代区间
%         - delta, epsilon分别为c,y的允许误差
%         - maxi为最大迭代次数
% 输出    - c为最终迭代结果
%         - err为c的误差值
%         - i为实际迭代次数
%         - yc为迭代结果c对应的y值
ya = f(a); yb = f(b);
if ya*yb > 0
    disp("Error: f(a)*f(b) > 0"); return;
end
c_mat = zeros(1,maxi);
for i = 1:maxi
    tmp = yb*(b-a)/(yb-ya);
    c = b - tmp;        % tmp = b-c，即bc之间的区间宽度
    c_mat(i) = c;
    a_c = c - a;
    yc = f(c);
    if yc == 0
        break;
    elseif ya*yc < 0
        b = c; yb = yc;
    else
        a = c; ya = yc;
    end
    if min(tmp,a_c) < delta || abs(yc) < epsilon
        break;
    end
end
err = min(tmp,a_c);
end