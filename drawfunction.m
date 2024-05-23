clc %清除命令窗口的内容
clear %清除工作空间的所有变量
clear all %清除工作空间的所有变量，函数，和MEX文件
clf %清除当前的Figure
close %关闭当前的Figure窗口
close all %关闭所有的Figure窗口

x=0:0.01:2; % 横坐标
y1=exp(x).*sin(10*pi*x); % 测试函数一
y2=exp(-(x-0.1).^2).*(sin(5.*pi.*x.^(3/4))).^6; % 测试函数二

figure
plot(x,y1,'-r');%画出函数图像
hold on;


[ymax,n]=max(y1);


plot(x(n),y1(n),'r*')%用plot画点，给特殊点做标记
title('测试函数一图像');
%加标注，text(x,y,str)中，x y是标注所在的位置，str是标注的内容
% a = 0.05:0.2:1.85
% b = exp(a).*sin(10*pi*a)
% for i=1:10
%     text(a(i)-0.1,b(i)+0.2,['P_m_a_x(',num2str(a(i)),',',num2str(b(i)),')'])
% end

text(x(n)-0.1,y1(n)+0.2,['P_m_a_x(',num2str(x(n)),',',num2str(y1(n)),')'])
xlabel('x');
ylabel('y');

figure
plot(x,y2,'-r');%画出函数图像
hold on;
[ymax,n]=max(y2);
plot(x(n),y2(n),'r*')%用plot画点，给特殊点做标记
title('测试函数二图像');
%加标注，text(x,y,str)中，x y是标注所在的位置，str是标注的内容
text(x(n)+0.02,y2(n)-0.04,['P_m_a_x(',num2str(x(n)),',',num2str(y2(n)),')'])
xlabel('x');
ylabel('y');
