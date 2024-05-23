function sigamincost() 
clc; 
%% �㷨��������
eranum=1000; % ��Ⱥ���� 500
popsize=200; % ��Ⱥ��ģ
pcross=0.8; %������ʣ�0.5-0.85 0.8
pmutation=0.18; %������ʣ�0.05-0.2 0.18

%% ����ϵ������
listnum = 10;%������Ŀ
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

PATH1 = [1 2 4 5 6 7 9 10]; % ����·��
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

%% ���1----��5-3���ϵ�������
omega
ri
alfai
betai
ai
bi

%% �㷨��ʼ
BestPop=zeros(eranum,listnum);%���ÿһ�������Ÿ���
Trace=zeros(eranum,1);%���ÿһ������С�ɱ�
% ��ʼ����Ⱥ
pop = zeros(popsize,listnum); 
pop = initpop(popsize,listnum,Cin,Ci0,Qi0,Ti0,Tin,Si0,omegai,Qc); 
% �������Ŵ�����
i = 1; 
while i <= eranum 
	for j=1:popsize 
		[realcost] = calcost(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4); %������ calcost
		value(j)=sum(realcost);%������Ӧ��,�����������ֵ����˺����������ֱ����Ϊ��Ӧ��ֵ
	end 
	[MinValue,Index]=min(value); 
	BestPop(i,:)=pop(Index,:);  %Ѱ����õ�
	Trace(i)=MinValue; 
	
	[selectpop]=selectchrom(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4); %�����������������������ѡ��
	[CrossOverPop]=CrossOver(selectpop,pcross,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);%����
	[NewPop]=Mutation(CrossOverPop,pmutation,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);%����
	pop=NewPop;%����
	i = i+1; 
    if(isval(pop,listnum,Qi0,Ci0,Cin,Si0,omegai,Qc)==0)
%         disp("�����Ч����")
    else
        disp("�����Ч����")
    end
end 

%% ���2----��5-4�Ż����
[Minfx,I]=min(Trace); 
Ti=BestPop(I,:); 
str1=sprintf('������ %d ��\n',I); 
disp(str1); 
disp('������ҵʱ��TiΪ\n'); 
Ti %������Ⱦɫ��
[Qi]=calQi(Ti)
[Ci]=calCi(Ti)
[Si]=calSi(Ti)

%% ����ɱ�(������ֵfx��
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

%% ����ؼ�·����Ҫ������
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

%% ���� Qi 
function Qi = calQi(realtime)
Qi = [0 0 0 0 0 0 0 0 0 0];
len = 10;
alfai = [0.3777    0.1156    0.1188    0.1458    0.3768    0.0071    0.1333    0.3697    0.1809    0.3761];
betai = [0.8298   -0.5189   -0.8471   -0.9275   -1.0493    2.0795   -1.5463   -5.0460   -1.2607   -3.2987];
for i =1:len 	
    Qi(i) = log(alfai(i)*realtime(i)+betai(i));
end
%% ���� Si 
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

%% ���� Ci 
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
%% ���������
%����һ�����߽�֮���������������߽�
function RN = CrtRandNum(LB,UB) 
scale = UB-LB; 
RN = round(rand()*scale); 
RN = RN+LB; 

%% ������Ӧ��ֵ
function fill = calfill(cost) 
fill = 1000/cost; 

%% ��ʼ����Ⱥ
function pop = initpop(popsize,listnum,Cin,Ci0,Qi0,Ti0,Tin,Si0,omegai,Qc); 
for i = 1:popsize 
	pop(i,:)= crtval(listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);
end

%% ������Ч����
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

%% �ж�һ�������Ƿ���Ч����Ч���� 1�����򷵻� 0
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

%% ѡ��
function [selectpop]=selectchrom(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4)
[m,n]=size(pop); 
for i=1:m 
	fit(i)=calfill(sum(calcost(pop,listnum,Cin,ri,Tin,PATH1,PATH2,PATH3,PATH4)));%�Ժ���ֵΪ��Ӧ��
end 
% ���㿹��Ũ��
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

% ������������
S = zeros(1,m); 
for i = 1:m 
	sm = 0; 
	for j = 1:m 
		sm = sm+((fit(i)-fit(j)) * (fit(i)-fit(j))); 
	end 
	S(i) = sqrt(sm);
end 
a = 0.5;%������������
P = zeros(1,m); 
for i= 1:m 
	P(i) = a*S(i)/sum(S)+(1-a)*exp(0-ro(i)); 
end 

selectprob=P/sum(P);%ѡ�����
prob=cumsum(selectprob);%�ۼ�ѡ�����
sumprob=[0 prob]; 
for i=1:m 
	selectpop(i,:)=pop(length(find(rand>=sumprob)),:); 
end 

%% �������
function dis = distance(point1,point2) 
len = length(point1); 
sum = 0; 
for i= 1:len 
	d = (point1(i) - point2(i))* (point1(i) - point2(i)); 
	sum = sum+d; 
end 
dis = sqrt(sum); 

%% ����
%�������--------------------------------------------------- 
%�Ը����еĸ��尴�ս�����н��н��泴�������Ը����н���
%�����ĸ�������������棬�����Ӵ����壬������н���ĸ�
%������Ϊ�����������Ϊż��
%--------------------------------------------------- 
function [NewPop]=CrossOver(OldPop,pcross,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc)%OldPop Ϊ������Ⱥ��pcross Ϊ�������
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
		[NewPop(y1(i+1),:),NewPop(y1(i+2),:)]=EqualCrossOver(OldPop(y1(i+1),:),OldPop(y1(i+2),:),listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc); 
	end 
end 
NewPop(y2,:)=OldPop(y2,:); 

%--------------------------------------------------- 
% ������Ⱦɫ����н������������Ϊ�������һ�� 2 ��Ⱦɫ�峤��
% ֮�������������Ⱦɫ���λ�õĻ�����н��������������Ⱦɫ
% ��Ϊ��Ч�ģ�����������Ⱦɫ��
%--------------------------------------------------- 
function [children1,children2]=EqualCrossOver(parent1,parent2,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc) 
L=length(parent1); 
hidecode=CrtRandNum(2,L); 
children1=zeros(1,L); 
children2=zeros(1,L); 
%����
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

%% ����
function [NewPop]=Mutation(OldPop,pmutation,listnum,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc) 
[m,n]=size(OldPop); 
r=rand(1,m); 
position=find(r<=pmutation); 
len=length(position); 
if len>=1
	for i=1:len 
		k=CrtRandNum(1,n); %���ñ����λ��
		OldPop(position(i),k) = CrtRandNum(Ti0(k),Tin(k)); 
		res= isval(OldPop(position(i),:),listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
		while res == 0 
			OldPop(position(i),:)= crtval(n,Qi0,Ti0,Tin,Ci0,Cin,Si0,omegai,Qc);
			res= isval(OldPop(position(i),:),listnum,Qi0,Ci0,Cin,Si0,omegai,Qc);
		end 
	end 
end 
NewPop=OldPop; % ���� 

