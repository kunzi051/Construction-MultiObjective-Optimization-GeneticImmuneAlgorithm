% basic.m 遗传算法
clear all; 
%function siga() 
clc; 
T1=clock; 
%--------------------------------------------------------------- 
LB = 0; % LB--------------自变量下限
UB = 2; % UB--------------自变量上限
eranum=100; % eranum----------种群的代数,取 100--1000 
popsize=50; % popsize---------每一代种群的规模；此可取 50--100 
pcross=0.8; % pcross-----------交叉的概率,此概率一般取 0.5--0.85 之间较好
pmutation=0.1; %pmutation------变异的概率,该概率一般取 0.05-0.2 左右较好
len = 20; %染色体长度，即二进制编码长度，它反应了数据的精度
%--------------------------------------------------------------- 
if find((LB-UB)>0) 
	error('数据输入错误,请重新输入(LB<UB):'); 
end 
s=sprintf('程序运行需要约%.4f秒钟时间,请稍等......',(eranum*popsize*40/(1000*50))); 
disp(s); 
% bounds=[LB;UB]';bits=[]; 
% precision=options(2);%由求解精度确定二进制编码长度
% bits=ceil(log2((bounds(:,2)-bounds(:,1))' ./ precision)); 
[Pop]=initpop(popsize,len);%初始化种群
[m,n]=size(Pop); 
pm0=pmutation; 
BestPop=zeros(eranum,n);%存放每一代的最优个体,二进制
Trace=zeros(eranum,2);%存放每一代的最大值，和最大值对应的变量值，实数
i=1; 
while i<=eranum 
	for j=1:m 
		temp = b2f(Pop(j,:),LB,UB); 
		value(j)=FUN(temp);%计算适应度,由于是求最大值，因此函数结果即可直接作为适应度值
	end 
	[MaxValue,Index]=max(value); 
	BestPop(i,:)=Pop(Index,:); 
	Trace(i,1)=MaxValue; 
	%--------------------------------- 
	%保存当前找到的最大值
	MaxValue = max(Trace(:,1)); 
% 	Trace(i,1)=(MaxValue+Trace(i,1))/2; 
	%--------------------------------- 
	Trace(i,2)=b2f(BestPop(i,:),LB,UB); 
	[selectpop]=selectchrom(Pop,LB,UB,MaxValue);%选择
	[CrossOverPop]=CrossOver(selectpop,pcross);%交叉
	[NewPop]=Mutation(CrossOverPop,pmutation);%变异
	Pop=NewPop;%更新
%	pmutation=pm0+(i^4)*(0.6-pm0)/(eranum^4); %随着种群向前进化，逐步增大变异率
%	p(i)=pmutation; 
	i=i+1; 
end 

%% 绘图
figure; 
t=1:eranum; 
plot(t,Trace(:,1)','-r'); % 线图，红色
title('测试函数二优化的遗传算法');xlabel('进化世代数(eranum)');ylabel('每一代最大值(maxfitness)'); 
[MaxFval,I]=max(Trace(:,1)); 
X=Trace(I,2); 
hold on; plot(I,MaxFval,'*r'); % 标记点
text(I+2,MaxFval,['F_m_a_x(',num2str(I),',',num2str(MaxFval),')']); % 标记坐标
str1=sprintf('进化到 %d 代 ,自变量为 %s 时,得本次求解的最优值 %f\n 对应染色体是：%s',I,num2str(X),MaxFval,num2str(BestPop(I,:))); 
disp(str1); 
disp('每代最优值对应的自变量为'); 
Trace(:,2) 
disp('每代的最优值是\n'); 
Trace(:,1) 
%figure(2);plot(t,p);%绘制变异值增大过程
T2=clock; 
CostTime=T2-T1; 
if CostTime(6)<0 
	CostTime(6)=CostTime(6)+60; CostTime(5)=CostTime(5)-1; 
end 
if CostTime(5)<0 
	CostTime(5)=CostTime(5)+60;CostTime(4)=CostTime(4)-1; 
end %像这种程序当然不考虑运行上小时啦
str2=sprintf(' 程序运行耗时 %d 小时 %d 分钟 %.4f 秒',CostTime(4),CostTime(5),CostTime(6)); 
disp(str2); 
%--------------------------------------------------- 
%初始化种群，采用二进制编码
%pop - 初始化后的种群
%popsize - 种群大小
%len - 每个染色体长度
function [pop]=initpop(popsize,len) 
pop(1,:)=zeros(1,len); 
for i=2:popsize-1 
	pop(i,:)=round(rand(1,len)); 
end
pop(popsize,:)=ones(1,len); 
end
%--------------------------------------------------- 
%解码,二进制到十进制的转换
% fval - 表征变量的十进制数
% LB - 变量取值下界
% UB - 变量取值上界
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
%% 选择操作
% selectpop -经过选择后的种群
% pop -原始种群
% LB -自变量下界
% UB -自变量上界
%--------------------------------------------------- 
function [selectpop]=selectchrom(pop,LB,UB,MaxValue)%采用轮盘赌选择方法
[m,n]=size(pop); 
for i=1:m 
	val(i) = b2f(pop(i,:),LB,UB);%求出各个体对应的十进制值
	res(i)=FUN ( val(i));%求出各个体对应的函数值
end 
delta = 1/((UB-LB)/2);%相似度抑制阈值
threshDis = 1/delta;%由抑制阈值计算出距离阈值
ro= zeros(1,m);
for i=1:m
	num = 0;
	for j= 1:m 
		dis = abs(val(i)-val(j));%由距离公式就为差的绝对值了
		if dis > threshDis %对小于相似度阈值的判断等价于对大于距离阈值的判断
			num=num+1; 
		end 
	end 
	ro(i) = num/m
end 
S = zeros(1,m);
for i =1:m 
	num = 0; 
	for j=1:m 
		num = num+abs(res(i) - res(j)); %抗体矢量距离计算公式
	end 
	S(i) = num; 
end 
alafa = 0.5; %常数调节因子，其值在【0，1】之间，可调节
SUM = sum(S); 
P=zeros(1,m); %存放各染色体的遗传概率
for i=1:m 
	P(i) = alafa*S(i)/SUM+(1-alafa)*exp(-ro(i)); 
end 
selectprob=P/sum(P);%选择概率
prob=cumsum(selectprob);%累计选择概率
sumprob=[0 prob]; 
for i=1:m 
	selectpop(i,:)=pop(length(find(rand>=sumprob)),:);
end 
end

%% 交叉变换
%输入变量：pop：二进制的父代种群数，pc：交叉的概率
%输出变量：newpop：交叉后的种群数
function [NewPop]=CrossOver(pop,pc)
[px,py] = size(pop);
NewPop = ones(size(pop));
for i = 1:2:px-1
    if(rand<pc)
        cpoint = round(rand*py);
        NewPop(i,:) = [pop(i,1:cpoint),pop(i+1,cpoint+1:py)];
        NewPop(i+1,:) = [pop(i+1,1:cpoint),pop(i,cpoint+1:py)];
    else
        NewPop(i,:) = pop(i,:);
        NewPop(i+1,:) = pop(i+1,:);
    end
end
end
%--------------------------------------------------- 
%% 变异操作
%--------------------------------------------------- 
function [NewPop]=Mutation(OldPop,pmutation) 
[m,n]=size(OldPop); 
r=rand(1,m); 
position=find(r<=pmutation); 
len=length(position); 
if len>=1 
	for i=1:len 
		k=unidrnd(n,1,1); %设置变异点数，一般设置 1 点，产生一个最大值为 n的服从统一分布的随机数
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
%单变量多峰的函数，需要测试哪个函数时，就将给函数命名为 FUN 
%--------------------------------------------------- 
function val = FUN(x) 
val = exp(-(x-0.1)^2)*(sin(5*pi*x^(3/4)))^6; 
% val = exp(x)*sin(10*pi*x); 
% val = x*sin(10*pi*x)+2; 
%--------------------------------------------------- 
end