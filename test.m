% test.m �Ż��������Ŵ��㷨
clear all; 
%function siga() 
clc; 
T1=clock; 
%--------------------------------------------------------------- 
LB = 0; % LB--------------�Ա�������
UB = 2; % UB--------------�Ա�������
eranum=100; % eranum----------��Ⱥ�Ĵ���,ȡ 100--1000 
popsize=50; % popsize---------ÿһ����Ⱥ�Ĺ�ģ���˿�ȡ 50--100 
pcross=0.8; % pcross-----------����ĸ���,�˸���һ��ȡ 0.5--0.85 ֮��Ϻ�
pmutation=0.1; %pmutation------����ĸ���,�ø���һ��ȡ 0.05-0.2 ���ҽϺ�
len = 20; %Ⱦɫ�峤�ȣ��������Ʊ��볤�ȣ�����Ӧ�����ݵľ���
%--------------------------------------------------------------- 
if find((LB-UB)>0) 
	error('�����������,����������(LB<UB):'); 
end 
s=sprintf('����������ҪԼ%.4f����ʱ��,���Ե�......',(eranum*popsize*40/(1000*50))); 
disp(s); 
% bounds=[LB;UB]';bits=[]; 
% precision=options(2);%����⾫��ȷ�������Ʊ��볤��
% bits=ceil(log2((bounds(:,2)-bounds(:,1))' ./ precision)); 
[Pop]=initpop(popsize,len);%��ʼ����Ⱥ
[m,n]=size(Pop); 
pm0=pmutation; 
BestPop=zeros(eranum,n);%���ÿһ�������Ÿ���,������
Trace=zeros(eranum,2);%���ÿһ�������ֵ�������ֵ��Ӧ�ı���ֵ��ʵ��
i=1; 
while i<=eranum 
	for j=1:m 
		temp = b2f(Pop(j,:),LB,UB); 
		value(j)=FUN(temp);%������Ӧ��,�����������ֵ����˺����������ֱ����Ϊ��Ӧ��ֵ
	end 
	[MaxValue,Index]=max(value); 
	BestPop(i,:)=Pop(Index,:); 
	Trace(i,1)=MaxValue; 
	%--------------------------------- 
	%���浱ǰ�ҵ������ֵ
% 	MaxValue = max(Trace(:,1)); 
% 	Trace(i,1)=MaxValue; 
	%--------------------------------- 
	Trace(i,2)=b2f(BestPop(i,:),LB,UB); 
	[selectpop]=selectchrom(Pop,LB,UB,MaxValue);%�����������������������ѡ��,�˴�����ڼ��Ŵ��㷨���˸Ľ�
	[CrossOverPop]=CrossOver(selectpop,pcross);%����
	[NewPop]=Mutation(CrossOverPop,pmutation);%����
	Pop=NewPop;%����
%	pmutation=pm0+(i^4)*(0.6-pm0)/(eranum^4); %������Ⱥ��ǰ�����������������
%	p(i)=pmutation; 
	i=i+1; 
end 

%% ��ͼ

figure; 
t=1:eranum; 
plot(t,Trace(:,1)','-r'); % ��ͼ����ɫ
title('���Ժ������Ż��������Ŵ��㷨');xlabel('����������(eranum)');ylabel('ÿһ�����ֵ(maxfitness)'); 
% ylim([2.5, 6.5])
[MaxFval,I]=max(Trace(:,1)); 
X=Trace(I,2); 
hold on; plot(I,MaxFval,'*r'); % ��ǵ�
text(I+2,MaxFval,['F_m_a_x(',num2str(I),',',num2str(MaxFval),')']); % �������
str1=sprintf('������ %d �� ,�Ա���Ϊ %s ʱ,�ñ�����������ֵ %f\n ��ӦȾɫ���ǣ�%s',I,num2str(X),MaxFval,num2str(BestPop(I,:))); 
disp(str1); 
disp('ÿ������ֵ��Ӧ���Ա���Ϊ'); 
Trace(:,2) 
disp('ÿ��������ֵ��\n'); 
Trace(:,1) 
T2=clock; 
CostTime=T2-T1; 
if CostTime(6)<0 
	CostTime(6)=CostTime(6)+60; CostTime(5)=CostTime(5)-1; 
end 
if CostTime(5)<0 
	CostTime(5)=CostTime(5)+60;CostTime(4)=CostTime(4)-1; 
end %
str2=sprintf(' �������к�ʱ %d Сʱ %d ���� %.4f ��',CostTime(4),CostTime(5),CostTime(6)); 
disp(str2); 
%--------------------------------------------------- 
%��ʼ����Ⱥ�����ö����Ʊ���
%pop - ��ʼ�������Ⱥ
%popsize - ��Ⱥ��С
%len - ÿ��Ⱦɫ�峤��
function [pop]=initpop(popsize,len) 
pop(1,:)=zeros(1,len); 
for i=2:popsize-1 
	pop(i,:)=round(rand(1,len)); 
end
pop(popsize,:)=ones(1,len); 
end
%--------------------------------------------------- 
%����,�����Ƶ�ʮ���Ƶ�ת��
% fval - ����������ʮ������
% LB - ����ȡֵ�½�
% UB - ����ȡֵ�Ͻ�
%--------------------------------------------------- 
function fval = b2f(bval,LB,UB) 
[X,Y] = size(bval); 
scale=(UB-LB)/(2^Y-1); 
sum = 0; 
for i=1:Y 
	sum = sum+bval(i)*2^(Y-i); 
end 
fval = LB+sum*scale; 
end
%--------------------------------------------------- 
%% ѡ�����
% selectpop -����ѡ������Ⱥ
% pop -ԭʼ��Ⱥ
% LB -�Ա����½�
% UB -�Ա����Ͻ�
%--------------------------------------------------- 
function [selectpop]=selectchrom(pop,LB,UB,MaxValue)%�������̶�ѡ�񷽷�
[m,n]=size(pop); 
for i=1:m 
	val(i) = b2f(pop(i,:),LB,UB);%����������Ӧ��ʮ����ֵ
	res(i)=FUN ( val(i));%����������Ӧ�ĺ���ֵ
end 
delta = 1/((UB-LB)/2);%���ƶ�������ֵ
threshDis = 1/delta;%��������ֵ�����������ֵ
ro= zeros(1,m);
for i=1:m
	num = 0;
	for j= 1:m 
		dis = abs(val(i)-val(j));%�ɾ��빫ʽ��Ϊ��ľ���ֵ��
		if dis > threshDis %��С�����ƶ���ֵ���жϵȼ��ڶԴ��ھ�����ֵ���ж�
			num=num+1; 
		end 
	end 
	ro(i) = num/m
end 
S = zeros(1,m);
for i =1:m 
	num = 0; 
	for j=1:m 
		num = num+abs(res(i) - res(j)); %����ʸ��������㹫ʽ
	end 
	S(i) = num; 
end 
alafa = 0.5; %�����������ӣ���ֵ�ڡ�0��1��֮�䣬�ɵ���
SUM = sum(S); 
P=zeros(1,m); %��Ÿ�Ⱦɫ����Ŵ�����
for i=1:m 
	P(i) = alafa*S(i)/SUM+(1-alafa)*exp(-ro(i)); 
end 
selectprob=P/sum(P);%ѡ�����
prob=cumsum(selectprob);%�ۼ�ѡ�����
sumprob=[0 prob]; 
for i=1:m 
	selectpop(i,:)=pop(length(find(rand>=sumprob)),:);
end 
end
%--------------------------------------------------- 
%% �������
%�Ը����еĸ��尴�ս�����н��н��泴�������Ը����н���
%�����ĸ�������������棬�����Ӵ����壬������н���ĸ�
%������Ϊ�����������Ϊż��
%--------------------------------------------------- 
function [NewPop]=CrossOver(OldPop,pcross)%OldPop Ϊ������Ⱥ��pcross Ϊ�������
[m,n]=size(OldPop); 
r=rand(1,m); 
y1=find(r<pcross);%�������ڽ�����к�
y2=find(r>=pcross); 
len=length(y1); 
if len>2&mod(len,2)==1%����������н����Ⱦɫ�������Ϊ�������������Ϊż��
	y2(length(y2)+1)=y1(len); 
	y1(len)=[]; 
end 
if length(y1)>=2 
	for i=0:2:length(y1)-2 
		[NewPop(y1(i+1),:),NewPop(y1(i+2),:)]=EqualCrossOver(OldPop(y1(i+1),:),OldPop(y1(i+2),:)); 
	end 
end 
NewPop(y2,:)=OldPop(y2,:); 
end
%--------------------------------------------------- 
%���þ��Ƚ��� ����
%�� 1�� 0 1 1 1 0 0 1 1 0 1 0 
%�� 2�� 1 0 1 0 1 1 0 0 1 0 1 
%���룺0 1 1 0 0 0 1 1 0 1 0
%������¸��壺
%�� 1�� 1 1 1 0 1 1 1 1 1 1 1 
%�� 2�� 0 0 1 1 0 0 0 0 0 0 0 
%--------------------------------------------------- 
function [children1,children2]=EqualCrossOver(parent1,parent2) 
L=length(parent1); 
hidecode=round(rand(1,L));%����������룬�� hidecode=[0 1 1 0 0 0 1 1 0 1 0]; 
children1=zeros(1,L); 
children2=zeros(1,L); 
children1(find(hidecode==1))=parent1(find(hidecode==1));%����Ϊ 1���� 1 Ϊ�� 1�ṩ����
children1(find(hidecode==0))=parent2(find(hidecode==0));%����Ϊ 0���� 2 Ϊ�� 1�ṩ����
children2(find(hidecode==1))=parent2(find(hidecode==1));%����Ϊ 1���� 2 Ϊ�� 2�ṩ����
children2(find(hidecode==0))=parent1(find(hidecode==0));%����Ϊ 0���� 1 Ϊ�� 2�ṩ����
end
%--------------------------------------------------- 
%% �������
%--------------------------------------------------- 
function [NewPop]=Mutation(OldPop,pmutation) 
[m,n]=size(OldPop); 
r=rand(1,m); 
position=find(r<=pmutation); 
len=length(position); 
if len>=1 
	for i=1:len 
		k=unidrnd(n,1,1); %���ñ��������һ������ 1 �㣬����һ�����ֵΪ n�ķ���ͳһ�ֲ��������
		for j=1:length(k) 
			if OldPop(position(i),k(j))==1 
				OldPop(position(i),k(j))=0; 
			else 
				OldPop(position(i),k(j))=1;
			end 
		end 
	end 
end 
NewPop=OldPop; 
end
%--------------------------------------------------- 
%--------------------------------------------------- 
%% ���������ĺ���
%��Ҫ�����ĸ�����ʱ���ͽ�����������Ϊ FUN 
%--------------------------------------------------- 
function val = FUN(x) 
val = exp(-(x-0.1)^2)*(sin(5*pi*x^(3/4)))^6; 
% val = exp(x)*sin(10*pi*x); 

%--------------------------------------------------- 
end