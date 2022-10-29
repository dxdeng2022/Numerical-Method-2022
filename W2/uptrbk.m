function X = uptrbk(A,B)
% 选主元的高斯消去方法(上三角变换和回代过程)
% 输入    - A是一个n*n非奇异矩阵
%         - B是一个n*1的向量
% 输出    - X是方程AX = B的解
N = size(A,1);
Aug = [A B];        % 增广矩阵
for p = 1:N-1
    % 选第p列的主元
    [~,j] = max(abs(Aug(p:N,p)));
    C = Aug(p,:);
    Aug(p,:) = Aug(p+j-1,:);
    Aug(p+j-1,:) = C;
    if Aug(p,p) == 0
        disp('Function uptrbk : A was singular matrix. No unique solution');
        break;
    end
    % 消去第p列对角线下的元素
    for k = p+1:N
        m = Aug(k,p)/Aug(p,p);
        Aug(k,p:N+1) = Aug(k,p:N+1)-m*Aug(p,p:N+1);
    end
end
X = backsub(Aug(1:N,1:N),Aug(1:N,N+1));        % 回代反算X
end


%% The function backsub(U,Y)
function X = backsub(U,Y)
% 输入    - U是一个n*n的上三角非奇异矩阵
%         - Y是n*1的列向量
% 输出    - X是线性方程UX=Y的解
N = length(Y);
X = zeros(N,1);
X(N) = Y(N)/U(N,N);
for k = N-1:-1:1
    X(k) = (Y(k)-U(k,k+1:N)*X(k+1:N))/U(k,k);
end
end
