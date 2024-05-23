clc %�������ڵ�����
clear %��������ռ�����б���
clear all %��������ռ�����б�������������MEX�ļ�
clf %�����ǰ��Figure
close %�رյ�ǰ��Figure����
close all %�ر����е�Figure����

x=0:0.01:2; % ������
y1=exp(x).*sin(10*pi*x); % ���Ժ���һ
y2=exp(-(x-0.1).^2).*(sin(5.*pi.*x.^(3/4))).^6; % ���Ժ�����

figure
plot(x,y1,'-r');%��������ͼ��
hold on;


[ymax,n]=max(y1);


plot(x(n),y1(n),'r*')%��plot���㣬������������
title('���Ժ���һͼ��');
%�ӱ�ע��text(x,y,str)�У�x y�Ǳ�ע���ڵ�λ�ã�str�Ǳ�ע������
% a = 0.05:0.2:1.85
% b = exp(a).*sin(10*pi*a)
% for i=1:10
%     text(a(i)-0.1,b(i)+0.2,['P_m_a_x(',num2str(a(i)),',',num2str(b(i)),')'])
% end

text(x(n)-0.1,y1(n)+0.2,['P_m_a_x(',num2str(x(n)),',',num2str(y1(n)),')'])
xlabel('x');
ylabel('y');

figure
plot(x,y2,'-r');%��������ͼ��
hold on;
[ymax,n]=max(y2);
plot(x(n),y2(n),'r*')%��plot���㣬������������
title('���Ժ�����ͼ��');
%�ӱ�ע��text(x,y,str)�У�x y�Ǳ�ע���ڵ�λ�ã�str�Ǳ�ע������
text(x(n)+0.02,y2(n)-0.04,['P_m_a_x(',num2str(x(n)),',',num2str(y2(n)),')'])
xlabel('x');
ylabel('y');
