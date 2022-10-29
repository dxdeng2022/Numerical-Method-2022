%% excrcise 1-1 -- 2022-09-24
% a = 0.2;
% for i = 1:2
%     a = a-0.1;
% end
% disp(a)


%% excrcise 1-10 -- 2022-09-24
% f = @(x) 4*x^3-2*x-6;
% df = @(x) 12*x^2-2;
% % 牛顿迭代法
% [p0,err,i,y,p_mat] = newton(f,df,10,5e-10,5e-10,1000);
% % 试值法
% [c,err2,i2,yc,c_mat] = regula(f,0,4,5e-10,5e-10,1000);
% % 绘图比较处理
% p_mat(i+2:i2) = p_mat(i+1);
% plot(p_mat(1:i2),'k','linewidth',1.5);
% hold on;
% plot(c_mat(1:i2),'r','linewidth',1.2);
% plot([0,i2],[p0,p0],'b--','linewidth',1.2);
% legend('Newton','False Position','approximation');
% set(gca,'fontname','Times New Roman','fontsize',14);
% grid on; hold off;


%% exercise 1-11 -- 2022-09-25
f = @(x) (x+4).^2.*(x+2).*(x-2).*(x-4).^3;
df = @(x) 2*(x+4)*(x+2)*(x-2)*(x-4)^3 + (x+4)^2*(2*x)*(x-4)^3 + 3*(x+4)^2*(x+2)*(x-2)*(x-4)^2;
p0 = 3;
[p0,err,i,yp0,p_mat] = newton(f,df,p0,5e-10,5e-10,1000);

x = -6:0.01:6;
y = f(x);
figure;
plot(x,y,'k','linewidth',1.5);
hold on; grid on;
yp = f(p_mat(1:i+1));
scatter(p_mat(1:i+1),yp,'ro','linewidth',1.5);
legend('Objective Function','Iteration Sequence');
set(gca,'fontname','Times New Roman','fontsize',14);
hold off;
