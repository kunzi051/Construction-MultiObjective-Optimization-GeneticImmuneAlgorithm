function sigamincost() 
clc; 
format long G;
%% 算法参数设置
eranum=2000; % 种群代数 1000;50
popsize=100; % 种群规模
pcross=0.5; %交叉概率，0.5-0.85 0.8
pmutation=0.18; %变异概率，0.05-0.2 0.18

%% 工程系数设置
listnum = 14;%工序数目
% global omegai Tin Ti0 Cin Ci0 Qi0 Sin Si0 epso1 epso2 epso3 epso4 k1 k2 k3 k4 c eta pai beta1 Tc Cc Cqc Csc Qc Sc
omegai = [0 0.22 0.24 0.22 0.2 0.35 0.61 0.59 0.62 0.58 0.54 0.01 0.15 0.46];
Tin = [5 28 35 30 10 90 12 30 15 2 40 21 22 16];
Ti0 = [4 25 32 28 9 45 11 28 13 1 36 20 20 15];
Cin = [251354     1693868     2340589     1860535     1350945     4413522      456628     1773846      653539      179605     4788317     1471842     3964281     6472140];
Ci0 = [254308     1703051     2349358     1866615     1353848     4248027      459676     1779714      659393      182589     4799941     1474801     3970329     6475113];
Qi0 = [ 0.8504    0.8635    0.8594    0.8665    0.8508    0.8750    0.8774    0.8589    0.8734    0.8621    0.8968    0.8538    0.8572    0.8511];
Sin = [0.9071    0.9537    0.9338    0.9378    0.9292    0.9334    0.9146    0.9297    0.9235    0.9134    0.9280    0.9257    0.9450    0.9094];
Si0 = [0.9016    0.9270    0.9122    0.9296    0.9255    0.9142    0.9100    0.9202    0.9143    0.9062    0.9166    0.9160    0.9278    0.9071];
epso1 = 0.3; epso2 = 0.3; epso3 = 0.25; epso4 = 0.15;
k1=1; k2=1; k3=0.5; k4=1; c=1.5; eta=0.03;
pai=1; beta1=13000; 
T=0;
Tc=250; % 工期
Cc=32000000; Cqc=16000000; Csc=960000; % 单位为元
Qc = 0.85; Sc = 0.9;

PATH1=[1 2 3 6 7 11 12]; % 工序路径，共12条
PATH2=[1 2 3 6 7 11 13];
PATH3=[1 2 3 6 8 11 12];
PATH4=[1 2 3 6 8 11 13];
PATH5=[1 2 3 6 9 10 11 12];
PATH6=[1 2 3 6 9 10 11 13];
PATH7=[1 2 4 5 6 7 11 12];
PATH8=[1 2 4 5 6 7 11 13];
PATH9=[1 2 4 5 6 8 11 12];
PATH10=[1 2 4 5 6 8 11 13];
PATH11=[1 2 4 5 6 9 10 11 12];
PATH12=[1 2 4 5 6 9 10 11 13];

for i=1:listnum

    ri(i) = (Ci0(i)-Cin(i))/((Tin(i)-Ti0(i))*(Tin(i)-Ti0(i)));
    alfai(i) = (exp(1)-exp(Qi0(i)))/(Tin(i)-Ti0(i));
    betai(i) = (exp(Qi0(i))*Tin(i)-exp(1)*Ti0(i))/(Tin(i)-Ti0(i));

    ai(i) = (exp(Sin(i))-exp(Si0(i)))/(Tin(i)-Ti0(i));
    bi(i) = (exp(Si0(i))*Tin(i)-exp(Sin(i))*Ti0(i))/(Tin(i)-Ti0(i));
    
end
omega = [0 0.0851    0.0768    0.0768    0.0729    0.1033    0.0812    0.0906    0.0713    0.0674    0.1127    0.0514    0.0552    0.0552]; 

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
pop = initpop(popsize,listnum,Cin,Ci0,Qi0,Ti0,Tin,Si0,omega,Qc); 
% 以下是遗传操作
i = 1; 
while i <= eranum 
	for j=1:popsize 
		[realcost] = calcost(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4,PATH5,PATH6,PATH7,PATH8,PATH9,PATH10,PATH11,PATH12); %主函数 calcost
		value(j)=sum(realcost);%计算适应度,由于是求最大值，因此函数结果即可直接作为适应度值
	end 
	[MinValue,Index]=min(value); 
	BestPop(i,:)=pop(Index,:);  %寻找最好的
	Trace(i)=MinValue; 
	
	[selectpop]=selectchrom(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4,PATH5,PATH6,PATH7,PATH8,PATH9,PATH10,PATH11,PATH12); %基于相似性适量距离的免疫选择
	[CrossOverPop]=CrossOver(selectpop,pcross,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc);%交叉
	[NewPop]=Mutation(CrossOverPop,pmutation,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc);%变异
	pop=NewPop;%更新
	i = i+1; 
    if(isval(pop,listnum,Qi0,Ci0,Cin,Si0,omega,Qc)==0)
        disp("输出无效个体")
    else
        disp("输出有效个体")
    end
end 

%% 输出2----表5-4优化结果
[Minfx,I]=min(Trace); 
Ti=BestPop(I,:); 
Ti(1)=5;
Ti(2)=28;
Ti(3)=35;
Ti(4)=30;
Ti(5)=10;
str1=sprintf('进化到 %d 代\n',I); 
disp(str1); 
disp('理想作业时间Ti为'); 
Ti %输出最佳染色体
[Qi]=calQi(Ti)
[Ci]=calCi(Ti)
[Si]=calSi(Ti)
[Si_sum] = calSi_sum(Si);
disp('最终安全水平为'); Si_sum(14)
Q = 0;
for i=1:14
    Q = Q+omega(i)*Qi(i);
end
disp('最终质量水平为'); Q

%% 计算成本(所求函数值fx）
function fx = calcost(Ti,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4,PATH5,PATH6,PATH7,PATH8,PATH9,PATH10,PATH11,PATH12)
C=0; Cq=0; Cs=0;
epso1 = 0.3; epso2 = 0.3; epso3 = 0.25; epso4 = 0.15;
k1=1; k2=1; k3=0.5; k4=1; c=1.5; eta=0.03;
pai=1; beta1=13000; 
T=0;
Tc=250; % 工期
Cc=32000000; Cqc=1600000; Csc=960000; % 单位为元

[Qi]=calQi(Ti);
[Ci]=calCi(Ti);
[Si]=calSi(Ti);
T = criticletime(Ti,PATH1,PATH2,PATH3,PATH4,PATH5,PATH6,PATH7,PATH8,PATH9,PATH10,PATH11,PATH12);

for i=1:listnum
    C=C+Cin(i)+ri(i)*(Ti(i)-Tin(i))*(Ti(i)-Tin(i))+beta1*T+pai*(T-Tc); % beta
    Cq=Cq+Ci(i)*((0.5*tan(pi/4)*Qi(i))^k1+0.5*(c*cot((pi/4)*(1+Qi(i)))^k2));
    Cs=Cs+eta*Ci(i)*(0.3*(tan(pi/2)*Si(i))^k3+0.7*(c*(1/Si(i))*(1/Si(i))-1)^k4);
end
fx=epso1*C/Cc+epso2*Cq/Cqc+epso3*Cs/Csc+epso4*T/Tc;

%% 计算关键路径需要的天数
function ct = criticletime(realtime,path1,path2,path3,path4,path5,path6,path7,path8,path9,path10,path11,path12) 
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
	sum4 =sum4+realtime(path4(i)); 
end 


len5 = length(path5); 
sum5 = 0; 
for i = 1:len5
	sum5 =sum5+realtime(path5(i)); 
end 

len6 = length(path6); 
sum6 = 0; 
for i = 1:len6 
	sum6 =sum6+realtime(path6(i)); 
end 

len7 = length(path7); 
sum7 = 0; 
for i = 1:len7 
	sum7 =sum7+realtime(path7(i)); 
end 

len8 = length(path8); 
sum8 = 0; 
for i = 1:len8 
	sum8 =sum8+realtime(path8(i)); 
end 

len9 = length(path9); 
sum9 = 0; 
for i = 1:len9 
	sum9 =sum9+realtime(path9(i)); 
end 

len10 = length(path10); 
sum10 = 0; 
for i = 1:len10 
	sum10 =sum10+realtime(path10(i)); 
end 

len11 = length(path11); 
sum11 = 0; 
for i = 1:len11 
	sum11 =sum11+realtime(path11(i)); 
end 
len12 = length(path12); 
sum12 = 0; 
for i = 1:len12 
	sum12 =sum12+realtime(path12(i)); 
end 
ct = max([sum1 sum2 sum3 sum4 sum5 sum6 sum7 sum8 sum9 sum10 sum11 sum12]); 

%% 计算 Qi 
function Qi = calQi(realtime)
Qi = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
len = 14;
alfai = [0.3777    0.1156    0.1188    0.1699    0.3768    0.0071    0.3136    0.1789    0.1616    0.3502    0.0666    0.3697    0.1809    0.3761];
betai = [0.8298   -0.5189   -1.4413   -2.3774   -1.0493    2.0795   -1.0454   -2.6475    0.2940    2.0180    0.0529   -5.0460   -1.2607   -3.2987];
Qi(1)=0.8504;
Qi(2)=0.8635;
Qi(3)=0.8594;
Qi(4)=0.8665;
Qi(5)=0.8508;
for i =6:len 	
    Qi(i) = log(alfai(i)*realtime(i)+betai(i));
end
%% 计算 Si 
function Si = calSi(realtime)
Si = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
len = 14;
ai = [0.0136    0.0228    0.0181    0.0104    0.0094    0.0011    0.0115    0.0120    0.0115    0.0179    0.0072    0.0244    0.0219    0.0057];
bi = [2.4092    1.9571    1.9099    2.2415    2.4390    2.4464    2.3583    2.1744    2.3451    2.4570    2.2427    2.0121    2.0902    2.3916];
Sin = [0.9071    0.9537    0.9338    0.9378    0.9292    0.9334    0.9146    0.9297    0.9235    0.9134    0.9280    0.9257    0.9450    0.9094];
Si(1) = 0.9071;
Si(2) = 0.9537;
Si(3) = 0.9338;
Si(4) = 0.9378;
Si(5) = 0.9292;
for i =6:len 	
    Si(i) = log(ai(i)*realtime(i)+bi(i));
end

%% 计算 Si_sum 含义：经过调整后的（公式3-13） ,Si_sum(14)为最终的安全水平
function Si_sum = calSi_sum(Si)
Si_sum = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
omegai = [0 0.22 0.24 0.22 0.2 0.35 0.61 0.59 0.62 0.58 0.54 0.01 0.15 0.46];
Sin = [0.9071    0.9537    0.9338    0.9378    0.9292    0.9334    0.9146    0.9297    0.9235    0.9134    0.9280    0.9257    0.9450    0.9094];
% 针对不同工序，m不同
for i =1:5 % ABCDE
	Si_sum(i) = (1-omegai(i)*(1-Sin(i)))*Si(i);
end

Si_sum(6) = (1-omegai(4)*(1-Sin(4))-omegai(5)*(1-Sin(5)))*Si(6); % F
for i =7:10 % GHIJ
	Si_sum(i) = (1-omegai(i)*(1-Sin(i)))*Si(i);
end
Si_sum(11) = (1-omegai(7)*(1-Sin(7))-omegai(8)*(1-Sin(8))-omegai(10)*(1-Sin(10)))*Si(11); % K
for i =12:13 % ML
	Si_sum(i) = (1-omegai(i)*(1-Sin(i)))*Si(i);
end
Si_sum(14) = (1-omegai(12)*(1-Sin(12))-omegai(13)*(1-Sin(13)))*Si(14); % N

%% 计算 Ci 
function Ci = calCi(realtime)
Ci = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
len = 14;
ri =  1.0e+03 * [2.9540    1.0203    0.9743    1.5200    2.9030   -0.0817    3.0480    1.4670    1.4635    2.9840    0.7265    2.9590    1.5120    2.9730];
Tin = [5 28 35 30 10 90 12 30 15 2 40 21 22 16];
Ti0 = [4 25 32 28 9 45 11 28 13 1 36 20 20 15];
Cin = [251354     1693868     2340589     1860535     1350945     4413522      456628     1773846      653539      179605     4788317     1471842     3964281     6472140];
Ci0 = [254308     1703051     2349358     1866615     1353848     4248027      459676     1779714      659393      182589     4799941     1474801     3970329     6475113];

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
function pop = initpop(popsize,listnum,Cin,Ci0,Qi0,Ti0,Tin,Si0,omega,Qc); 
for i = 1:popsize 
	pop(i,:)= crtval(listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc);
end

%% 产生有效个体
function person = crtval(listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc) 
person(1) = 5;
person(2) = 28;
person(3) = 35;
person(4) = 30;
person(5) = 10;
for i =6:listnum 
	person(i) = CrtRandNum(Ti0(i),Tin(i)); 
end
res = isval(person,listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
while res==0
    person(1) = 5;
    person(2) = 28;
    person(3) = 35;
    person(4) = 30;
    person(5) = 10;
	for i =6:listnum 
		person(i) = CrtRandNum(Ti0(i),Tin(i)); 
	end 
	res = isval(person,listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
end

%% 判断一个个体是否有效，有效返回 1，否则返回 0
function res = isval(person,listnum,Qi0,Ci0,Cin,Si0,omega,Qc)
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
    Q = Q+omega(i)*Qi(i);
end
if Q<Qc
    res = 0;
end

%% 选择
function [selectpop]=selectchrom(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4,PATH5,PATH6,PATH7,PATH8,PATH9,PATH10,PATH11,PATH12)
[m,n]=size(pop); 
for i=1:m 
	fit(i)=calfill(sum(calcost(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4,PATH5,PATH6,PATH7,PATH8,PATH9,PATH10,PATH11,PATH12)));%以函数值为适应度
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
function [NewPop]=CrossOver(OldPop,pcross,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc)%OldPop 为父代种群，pcross 为交叉概率
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
		[NewPop(y1(i+1),:),NewPop(y1(i+2),:)]=EqualCrossOver(OldPop(y1(i+1),:),OldPop(y1(i+2),:),listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc); 
	end 
end 
NewPop(y2,:)=OldPop(y2,:); 

%--------------------------------------------------- 
% 对两条染色体进行交叉操作，具体为随机产生一个 2 到染色体长度
% 之间的数，对两条染色体该位置的基因进行交换，如果产生的染色
% 体为无效的，则重新生成染色体
%--------------------------------------------------- 
function [children1,children2]=EqualCrossOver(parent1,parent2,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc) 
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
res= isval(children1,listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
while res == 0 
	children1= crtval(L,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc);
	res= isval(children1,listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
end 
res= isval(children2,listnum,Qi0,Ci0,Cin,Si0,omega,Qc); 
while res == 0 
	children2= crtval(L,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc);
	res= isval(children2,listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
end 

%% 变异
function [NewPop]=Mutation(OldPop,pmutation,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc) 
[m,n]=size(OldPop); 
r=rand(1,m); 
position=find(r<=pmutation); 
len=length(position); 
if len>=1
	for i=1:len 
		k=CrtRandNum(1,n); %设置变异点位置
		OldPop(position(i),k) = CrtRandNum(Ti0(k),Tin(k)); 
		res= isval(OldPop(position(i),:),listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
		while res == 0 
			OldPop(position(i),:)= crtval(n,Qi0,Ti0,Tin,Ci0,Cin,Si0,omega,Qc);
			res= isval(OldPop(position(i),:),listnum,Qi0,Ci0,Cin,Si0,omega,Qc);
		end 
	end 
end 
NewPop=OldPop; % 更新 

