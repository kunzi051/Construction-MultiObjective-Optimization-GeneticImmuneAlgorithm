function sigamincost() 
clc; 
%% 算法参数设置
eranum=1000; % 种群代数 500
popsize=200; % 种群规模
pcross=0.8; %交叉概率，0.5-0.85 0.8
pmutation=0.18; %变异概率，0.05-0.2 0.18

%% 工程系数设置
listnum = 10;%工序数目
omegai = [0 0.22 0.24 0.22 0.2 0.35 0.54 0.01 0.15 0.46];
Tin = [5 28 30 25 10 90 32 21 22 16];
Ti0 = [4 25 27 23 9 45 30 20 20 15];
Cin = [125677 557495 938964 696408 387273 1915355 2144194 478760 1715355 3001871];
Ci0 = [130838 573773 960342 707344 392403 1715355 2157488 484314 1725359 3004054];
Qi0 = [0.8504 0.8635 0.8594 0.8865 0.8508 0.8750 0.8968 0.8538 0.8572 0.8511];
Sin = [0.9071 0.9537 0.9338 0.9378 0.9292 0.9134 0.9280 0.9257 0.945 0.9094];
Si0 = [0.9016 0.9270 0.9122 0.9296 0.9255 0.9042 0.9166 0.9160 0.9278 0.9071];
epso1 = 0.3; epso2 = 0.3; epso3 = 0.25; epso4 = 0.15;
k1=1; k2=1; k3=0.5; k4=1; c=1.5; eta=0.03;
pai=1; beta1=13000; T=0;
Tc=216; Cc=1200; Cqc=60; Csc=36;
Qc = 0.85; Sc = 0.9;

PATH1 = [1 2 4 5 6 7 9 10]; % 工序路径
PATH2 = [1 2 3 6 7 9 10];
PATH3 = [1 2 3 6 7 8 10];
PATH4 = [1 2 4 5 6 7 8 10]; 

for i=1:listnum

    ri(i) = (Ci0(i)-Cin(i))/((Tin(i)-Ti0(i))*(Tin(i)-Ti0(i)));
    alfai(i) = (exp(1)-exp(Qi0(i)))/(Tin(i)-Ti0(i));
    betai(i) = (exp(Qi0(i))*Tin(i)-exp(1)*Ti0(i))/(Tin(i)-Ti0(i));

    ai(i) = (exp(Sin(i))-exp(Si0(i)))/(Tin(i)-Ti0(i));
    bi(i) = (exp(Si0(i))*Tin(i)-exp(Sin(i))*Ti0(i))/(Tin(i)-Ti0(i));
    
end
omega = [0 77 69.5 69.5 66 93.5 102 46.5 50 50]/(0+77+69.5+69.5+66+93.5+102+46.5+50+50);

%% 输出1----表5-3相关系数计算表
omega
ri
alfai
betai
ai
bi

%% 算法开始
BestPop=zeros(eranum,listnum);%存放每一代的最优个体
Trace=zeros(eranum,1);%存放每一代的最小成本
% 初始化种群
pop = zeros(popsize,listnum); 
pop = initpop(popsize,listnum,Cin,Ci0,Qi0,Ti0,Tin,Si0,omegai,Qc); 
% 以下是遗传操作
i = 1; 
while i <= eranum 
	for j=1:popsize 
		[realcost] = calcost(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4); %主函数 calcost
		value(j)=sum(realcost);%计算适应度,由于是求最大值，因此函数结果即可直接作为适应度值
	end 
	[MinValue,Index]=min(value); 
	BestPop(i,:)=pop(Index,:);  %寻找最好的
	Trace(i)=MinValue; 
	
	[selectpop]=selectchrom(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4); %基于相似性适量距离的免疫选择
	[CrossOverPop]=CrossOver(selectpop,pcross,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);%交叉
	[NewPop]=Mutation(CrossOverPop,pmutation,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);%变异
	pop=NewPop;%更新
	i = i+1; 
    if(isval(pop,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc)==0)
%         disp("输出无效个体")
    else
        disp("输出有效个体")
    end
end 

%% 输出2----表5-4优化结果
[Minfx,I]=min(Trace); 
Ti=BestPop(I,:); 
str1=sprintf('进化到 %d 代\n',I); 
disp(str1); 
disp('理想作业时间Ti为\n'); 
Ti %输出最佳染色体
[Qi]=calQi(Ti)
[Ci]=calCi(Ti)
[Si]=calSi(Ti)

%% 计算成本(所求函数值fx）
function fx = calcost(Ti,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4)
C=0; Cq=0; Cs=0;
epso1 = 0.3; epso2 = 0.3; epso3 = 0.25; epso4 = 0.15;
k1=1; k2=1; k3=0.5; k4=1; c=1.5; eta=0.03;
pai=1; beta1=13000; T=0;
Tc=216; Cc=1200; Cqc=60; Csc=36;
[Qi]=calQi(Ti);
[Ci]=calCi(Ti);
[Si]=calSi(Ti);
T = criticletime(Ti,PATH1,PATH2,PATH3,PATH4);
% for i=1:listnum
%     T=T+Ti(i);
% end
for i=1:listnum
    C=C+Cin(i)+ri(i)*(Ti(i)-Tin(i))*(Ti(i)-Tin(i))+beta1*T+pai*(T-Tc); % beta
    Cq=Cq+Ci(i)*((0.5*tan(pi/4)*Qi(i))^k1+0.5*(c*cot((pi/4)*(1+Qi(i)))^k2));
    Cs=Cs+eta*Ci(i)*(0.3*(tan(pi/2)*Si(i))^k3+0.7*(c*(1/Si(i))*(1/Si(i))-1)^k4);
end
fx=epso1*C/Cc+epso2*Cq/Cqc+epso3*Cs/Csc+epso4*T/Tc;

%% 计算关键路径需要的天数
function ct = criticletime(realtime,path1,path2,path3,path4) 
len1 = length(path1); 
sum1 = 0; 
for i = 1:len1 
	sum1 =sum1+realtime(path1(i)); 
end 
len2 = length(path2); 
sum2 = 0; 

for i = 1:len2 
	sum2 =sum2+realtime(path2(i)); 
end 
len3 = length(path3); 
sum3 = 0; 
for i = 1:len3 
	sum3 =sum3+realtime(path3(i)); 
end 
len4 = length(path4); 
sum4 = 0; 
for i = 1:len4 
	sum3 =sum3+realtime(path4(i)); 
end 
ct = min([sum1 sum2 sum3 sum4]); 

%% 计算 Qi 
function Qi = calQi(realtime)
Qi = [0 0 0 0 0 0 0 0 0 0];
len = 10;
alfai = [0.3777    0.1156    0.1188    0.1458    0.3768    0.0071    0.1333    0.3697    0.1809    0.3761];
betai = [0.8298   -0.5189   -0.8471   -0.9275   -1.0493    2.0795   -1.5463   -5.0460   -1.2607   -3.2987];
for i =1:len 	
    Qi(i) = log(alfai(i)*realtime(i)+betai(i));
end
%% 计算 Si 
function Si = calSi(realtime)
Si1 = [0 0 0 0 0 0 0 0 0 0];
Si = [0 0 0 0 0 0 0 0 0 0];
omegai = [0    0.1234    0.1114    0.1114    0.1058    0.1498    0.1635    0.0745    0.0801    0.0801];
len = 10;
ai = [0.0136    0.0228    0.0181    0.0104    0.0094    0.0005    0.0143    0.0244    0.0219    0.0057];
bi = [2.4092    1.9571    2.0005    2.2936    2.4390    2.4471    2.0707    2.0121    2.0902    2.3916];
Sin =[0.9071 0.9537 0.9338 0.9378 0.9292 0.9134 0.9280 0.9257 0.945 0.9094];
for i =1:len 	
    Si1(i) = log(ai(i)*realtime(i)+bi(i));
end
for i =1:8
	Si(i) = (1-omegai(i)*(1-Sin(i)))*Si1(i);
end
Si(9) = (1-omegai(9)*(1-Sin(9))-omegai(8)*(1-Sin(8)))*Si1(9);
Si(10) = (1-omegai(10)*(1-Sin(10)))*Si1(10);

%% 计算 Ci 
function Ci = calCi(realtime)
Ci = [0 0 0 0 0 0 0 0 0 0];
len = 10;
ri =  1.0e+03 * [5.1610    1.8087    2.3753    2.7340    5.1300   -0.0988    3.3235    5.5540    2.5010    2.1830];
Tin = [5 28 30 25 10 90 32 21 22 16];
Ti0 = [4 25 27 23 9 45 30 20 20 15];
Cin = [125677 557495 938964 696408 387273 1915355 2144194 478760 1715355 3001871];
Ci0 = [130838 573773 960342 707344 392403 1715355 2157488 484314 1725359 3004054];
for i =1:len 	
    Ci(i) = Cin(i)+ri(i)*(realtime(i)-Tin(i))*(realtime(i)-Tin(i));  
end 
%% 生成随机数
%产生一个两边界之间的随机数，包括边界
function RN = CrtRandNum(LB,UB) 
scale = UB-LB; 
RN = round(rand()*scale); 
RN = RN+LB; 

%% 计算适应度值
function fill = calfill(cost) 
fill = 1000/cost; 

%% 初始化种群
function pop = initpop(popsize,listnum,Cin,Ci0,Qi0,Ti0,Tin,Si0,omegai,Qc); 
for i = 1:popsize 
	pop(i,:)= crtval(listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);
end

%% 产生有效个体
function person = crtval(listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc) 
for i =1:listnum 
	person(i) = CrtRandNum(Ti0(i),Tin(i)); 
end
res = isval(person,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
while res==0 
	for i =1:listnum 
		person(i) = CrtRandNum(Ti0(i),Tin(i)); 
	end 
	res = isval(person,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
end

%% 判断一个个体是否有效，有效返回 1，否则返回 0
function res = isval(person,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc)
res = 1;
[Qi] = calQi(person); 
[Si] = calSi(person); 
[Ci] = calCi(person); 

if Qi<Qi0
    res = 0;
end
if Qi>=1
    res = 0;
end
if Qi<0.85
    res = 0;
end
if Ci<Cin
    res = 0;
end
if Ci>Ci0
    res = 0;
end
if Si<Si0
    res = 0;
end
if Si<0.9
    res = 0;
end
if Si>=1
    res = 0;
end
Q = 0;
for i=1:listnum
    Q = Q+omegai(i)*Qi(i);
end
if Q<Qc
    res = 0;
end

%% 选择
function [selectpop]=selectchrom(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4)
[m,n]=size(pop); 
for i=1:m 
	fit(i)=calfill(sum(calcost(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4)));%以函数值为适应度
end 
% 计算抗体浓度
ro = zeros(1,m); 
for i = 1:m 
	r = 0; 
	for j = 1:m 
		dis1 = distance(pop(i,:),pop(j,:)); 
		dis2 = distance(pop(1,:),pop(j,:)); 
		if dis1>dis2 
			r = r+1; 
		end 
	end 
	ro(i) = r/n; 
end 

% 计算适量距离
S = zeros(1,m); 
for i = 1:m 
	sm = 0; 
	for j = 1:m 
		sm = sm+((fit(i)-fit(j)) * (fit(i)-fit(j))); 
	end 
	S(i) = sqrt(sm);
end 
a = 0.5;%常数调节因子
P = zeros(1,m); 
for i= 1:m 
	P(i) = a*S(i)/sum(S)+(1-a)*exp(0-ro(i)); 
end 

selectprob=P/sum(P);%选择概率
prob=cumsum(selectprob);%累计选择概率
sumprob=[0 prob]; 
for i=1:m 
	selectpop(i,:)=pop(length(find(rand>=sumprob)),:); 
end 

%% 两点距离
function dis = distance(point1,point2) 
len = length(point1); 
sum = 0; 
for i= 1:len 
	d = (point1(i) - point2(i))* (point1(i) - point2(i)); 
	sum = sum+d; 
end 
dis = sqrt(sum); 

%% 交叉
%交叉操作--------------------------------------------------- 
%对父代中的个体按照交叉比列进行交叉炒作，即对父代中交叉
%比例的个体进行两两交叉，产生子代个体，如果进行交叉的父
%代个数为奇数，则调整为偶数
%--------------------------------------------------- 
function [NewPop]=CrossOver(OldPop,pcross,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc)%OldPop 为父代种群，pcross 为交叉概率
[m,n]=size(OldPop);
r=rand(1,m); 
y1=find(r<pcross);%返回用于交叉的行号
y2=find(r>=pcross); 
len=length(y1); 
if len>2&mod(len,2)==1%如果用来进行交叉的染色体的条数为奇数，将其调整为偶数
	y2(length(y2)+1)=y1(len); 
	y1(len)=[]; 
end 
if length(y1)>=2 
	for i=0:2:length(y1)-2 
		[NewPop(y1(i+1),:),NewPop(y1(i+2),:)]=EqualCrossOver(OldPop(y1(i+1),:),OldPop(y1(i+2),:),listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc); 
	end 
end 
NewPop(y2,:)=OldPop(y2,:); 

%--------------------------------------------------- 
% 对两条染色体进行交叉操作，具体为随机产生一个 2 到染色体长度
% 之间的数，对两条染色体该位置的基因进行交换，如果产生的染色
% 体为无效的，则重新生成染色体
%--------------------------------------------------- 
function [children1,children2]=EqualCrossOver(parent1,parent2,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc) 
L=length(parent1); 
hidecode=CrtRandNum(2,L); 
children1=zeros(1,L); 
children2=zeros(1,L); 
%交叉
for i = 1:hidecode-1 
	children1(i) = parent1(i); 
	children2(i) = parent2(i);
end 
for i = hidecode:L 
	children1(i) = parent2(i); 
	children2(i) = parent1(i); 
end 
res= isval(children1,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
while res == 0 
	children1= crtval(L,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);
	res= isval(children1,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
end 
res= isval(children2,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc); 
while res == 0 
	children2= crtval(L,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);
	res= isval(children2,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
end 

%% 变异
function [NewPop]=Mutation(OldPop,pmutation,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc) 
[m,n]=size(OldPop); 
r=rand(1,m); 
position=find(r<=pmutation); 
len=length(position); 
if len>=1
	for i=1:len 
		k=CrtRandNum(1,n); %设置变异点位置
		OldPop(position(i),k) = CrtRandNum(Ti0(k),Tin(k)); 
		res= isval(OldPop(position(i),:),listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
		while res == 0 
			OldPop(position(i),:)= crtval(n,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);
			res= isval(OldPop(position(i),:),listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
		end 
	end 
end 
NewPop=OldPop; % 更新 

