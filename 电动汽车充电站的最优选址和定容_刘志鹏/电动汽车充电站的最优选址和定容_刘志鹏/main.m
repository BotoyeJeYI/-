clear ; 
clc; 
close all
%% 基础数据
%    1    2     3         4        5
%b=[需求坐标，需求负荷，成本系数，区域]
%充电需求点坐标
b(:,1) =1.0e+003.*([0.1092, 0.1197,0.2578,0.4259, 0.1257,0.2803,0.4439,0.5505,0.5610,0.5700,0.1332,0.3013,0.4559,0.5850,0.1452,0.3163,0.4739,0.5880,0.3193,0.4784,0.6015,0.6736,0.8657,1.0308,0.6781,0.8327,1.0188,0.6811,0.8192,1.0083,0.6691,0.6976,0.8191,1.0098]');
b (:,2)=1.0e+003.*[0.1348;0.2399;0.1724;0.1739;0.3420;0.3375;0.3360;0.1108;0.2024;0.3075;0.4591;0.4455;0.4380;0.4260;0.6092;0.5341;0.5341;0.5341;0.6407;0.6452; 0.6467;0.1574 ;0.1649;0.1634;0.2924;0.2909;0.2939;0.4425; 0.4576;0.4546;0.5431;0.6392;0.6377;0.6347];
%充电需求点常规电力负荷点负荷b(:,3)（kW）
b(:,3)=[2480;2480;8680;11400;890;2340;4160;560;1670;5010;2670;8280;7400;1430;7500;4840;3400;4290;3840;3680;2560;7000;14800;8960;3160;7000;5000;2280;10360;10000;760;6000;7040;5600];

%集中充电站坐标,两个，各区域一个
bcs=[ 937.7296  379.5010;
  310.3141  238.4076];
%选址数量
Tn=7;   

%% 成本参数
na=4500;    %电动汽车数量
alp=0.1;  
b(:,4)=round(alp.*b(:,3)./sum(b(:,3)).*na); %每个需求点平均负荷
b(23,4)=37;  
ns=4;
mui=0.6;   
Nchz=round(mui.*sum(b(:,4))./ns);
bm=1.0e+003*[0.0086,0.0088;1.1734,0.0088;1.1734,0.7412;0.0086,0.7412;0.0086,0.0088];
BL=sqrt(8.2*1.0e6./((max(bm(:,1))-min(bm(:,1)))*(max(bm(:,2))-min(bm(:,2)))));%BL为图坐标与实际坐标的比例，为固定参数

%% 区域划分
Area1=1.0e+003 *[0.9377   -1.0860;
1.1734    0.0088;
1.1734    0.7412;
0.3103    1.7040;
0.9377   -1.0860];
Area1=[Area1,1.*ones(size(Area1,1),1)]; %其中size(Area1,1)输出的是Area1的行数
Area2=1.0e+003 *[0.0086    0.0088;
0.9377   -1.0860;
0.3103    1.7040;
0.0086    0.7412;
0.0086    0.0088];
Area2=[Area2,2.*ones(size(Area2,1),1)];
%区域分界
vv=[Area1;Area2];   %10*3数组
for k=1:size(bcs,1) %size(bcs,1)输出bcs的行数即2,即k=1:2
    Ai=find(vv(:,3)==k);    %在vv的第三列查找等于k的元素返回索引值,即判别vv元素属于区域几，并作为索引暂时存入Ai中
    xx=vv(Ai,1);            %横坐标，充电需求点负荷，等于k的点 列向量
    yy=vv(Ai,2);            %纵坐标
    kk=convhull(xx,yy);     %计算凸包，kk是一个列向量
    %in = inpolygon(x,y,xv,yv)%注意xv,yv构成了多边形边界。x,y对应的是单点坐标，判断是否在多边形内。
    %返回结果为逻辑logical类型（不是数字类型哦），如果在对应的就返回1，否则为0。
    in=inpolygon(b(:,1),b(:,2),xx(kk),yy(kk));
    b(in,5)=k; 
end
%区域归属。1：成本系数之和，2：mui*成本系数之和/ns的值取整，3：区域
Ep=[];
for i=1:size(bcs,1)     %2个区域，即i=1：2
    gb=b(b(:,5)==i,:);  %把区域内的数据暂存到gb
    Ep=[Ep;[sum(gb(:,4)),round(mui.*sum(gb(:,4))./ns),i]];  %区域归属的数据计算及储存
end

%% 粒子群算法参数

PopSize=20; %种群数量
MaxIter=300;  %最大迭代次数
%学习因子
c1s=2.5;    %个体学习速率最大值
c2s=0.5;    %社会学习速率最大值
c1e=0.5;    %个体学习速率最小值
c2e=2.5;    %社会学习速率最小值
w_start=0.9;   %权值最大值
w_end=0.4;     %权值最小值
Iter=1;        

%上下限约束（即划定区域）
xmax=max(Area1(:,1));   %区域1最大长
xmin=min(Area1(:,1));   %区域1最小长
ymax=max(Area1(:,2));   %区域2最大宽
ymin=min(Area1(:,2));   %区域2最小宽
x = xmin + (xmax-xmin).*rand(Tn,PopSize);   %生成x方向的粒子，Tn个选址+PopSize个种群
y = ymin + (ymax-ymin).*rand(Tn,PopSize);   %生成y方向的粒子
%% 初始化
X=[x;y];    %合并xy粒子，形成一个完整的个体,前7行表示X坐标，后7行表示Y坐标，这样以列向量表示，供函数使用
V=rand(Tn*2,PopSize);   %随机生成粒子速度
Vmax=0.1*max((xmax-xmin),(ymax-ymin));      %最大粒子速度
inAr1=find(b(:,5)==1);                      %找到区域1的编号（索引）并存储到inAr1中
bb=[b(inAr1,1:2),b(inAr1,4)];               %bb暂存需求坐标和成本系数
%初始寻优
for pk=1:1:PopSize
    [FX(pk),~,~,~,~,~,~,~,~,~]=VorCostCDEV(X(1:Tn,pk),X(Tn+1:end,pk),bb,bcs(1,:),BL);   %求适应度值,~代表不输出某个参数，所以只输出第一个参数。X(1:Tn,pk)表示pk个种群的Tn个最优选址的横坐标,X(Tn+1:end,pk)表示pk个种群的Tn个最优选址的纵坐标
end

%输出初始最优
PBest=X;            %PBest初始最优坐标
FPBest=FX;          %FPBest初始最优目标函数值


[Fgbest,r]=min(FX); %最优目标函数值的初始最小值【最小值，第几列（即第几个种群）】
Fm(Iter)=Fgbest;    %把最优目标函数值的初始最小值存入Fm中
CF=Fgbest;          %把最优目标函数值的初始最小值存入CF中
Best=X(:,r);        %初始最优种群的坐标存入Best（列向量，前7横坐标，后7纵坐标）
FBest(Iter)=Fgbest; %把最优目标函数值的初始最小值存入FBest（1）中

FgNum=0;
%% 粒子群寻优
while (Iter<=MaxIter)   %粒子群算法
    Iter=Iter+1;        %迭代次数
    w_now=((w_start-w_end)*(MaxIter-Iter)/MaxIter)+w_end;   %惯性权重 
    A=repmat(X(:,r),1,PopSize);     
    R1=rand(Tn*2,PopSize);
    R2=rand(Tn*2,PopSize);   
    c1=c1e+(c1s-c1e)*(1-acos(-2*Iter/(MaxIter+1)+1)/pi);
    c2=c2e+(c2s-c2e)*(1-acos(-2*Iter/(MaxIter+1)+1)/pi);    
    V=w_now*V+c1*R1.*(PBest-X)+c2*R2.*(A-X);                %粒子速度更新公式
    changeRows=V>Vmax;
    V(changeRows)=Vmax;   
    changeRows=V<-Vmax;
    V(changeRows)=-Vmax;    
    X=X+1.0*V;  
   for pk=1:1:PopSize
       [FX(pk),~,~,~,~,~,~,~,~,~]=VorCostCDEV(X(1:Tn,pk),X(Tn+1:end,pk),bb,bcs(1,:),BL);    %求适应度
   end    
    P=FX<FPBest;    
    FPBest(P)=FX(P);  
    PBest(:,P)=X(:,P);     
    [Fgbest,r]=min(FPBest);
    Fm(Iter)=Fgbest;        
    
    if Fgbest<CF 
        [FBest,gr]=min(FPBest);  
        Best=PBest(:,gr);    
        CF=Fgbest;  
        FgNum=0;
    else
        FgNum=FgNum+1; 
    end
    
    
    if FgNum>10    
        SubX=r;
        while SubX==r || SubX==0
            SubX=round(rand*(PopSize));
        end
        r=SubX;
        FgNum=0;
    end
    
end 

%% 输出最优值
 %坐标
x =Best(1:Tn);
y =Best(Tn+1:end);
[f1,f2,f3,f4,Num_chI,num_Epc,f5,f6,f7,f8]=VorCostCDEV(x,y,bb,bcs(1,:),BL)
FBest  
Best 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 作图
%% 迭代图
figure(1)
Iter_t=1:1:MaxIter+1;
plot(Iter_t,Fm,'.-')
legend('每次迭代值')
title('迭代图')

%% 选址图
figure(2)
a=imread('T3.Png');
imshow(a);   
hold on
[vxT,vyT]=VoronoiT(bcs(:,1),bcs(:,2),0);  
%集中式充电站
plot(bcs(:,1),bcs(:,2),'gs','linewidth',12);
%区域分界
plot(vxT,vyT,'r -','linewidth',3);   
%需求点
plot(b(:,1),b(:,2),'k*','linewidth',5)  ;
%周边
plot(bm(:,1),bm(:,2),'k-','linewidth',3) 
%需求点排序
for k=1:length(b(:,4))
    str=num2str(b(k,4));
    text(b(k,1)-15,b(k,2)+25,str,'FontSize',23,'color','black');
end
axis equal
 [vx,vy]=voronoi(x,y);
 plot(x,y,'b^','linewidth',3); 
 %最优充电站排序
for k=1:length(x)
    str=num2str(k);
    text(x(k),y(k),str,'FontSize',20,'color','red');
 end
axis([min(bm(:,1))-3 max(bm(:,1))+3 min(bm(:,2))-3 max(bm(:,2))+3])
legend('集中充电站','区域分界','充电需求点','周边界限','充电站最优选址区')
title('规划选址图')

%% 选址图
figure(3)
a=imread('T3.bmp');
imshow(a);   
hold on
[vxT,vyT]=VoronoiT(bcs(:,1),bcs(:,2),0);  
%集中式充电站
plot(bcs(:,1),bcs(:,2),'gs','linewidth',12);
%区域分界
plot(vxT,vyT,'r -','linewidth',3);   
%需求点
plot(b(:,1),b(:,2),'k*','linewidth',5)  ;
%周边
plot(bm(:,1),bm(:,2),'k-','linewidth',3) 
%需求点排序
for k=1:length(b(:,4))
    str=num2str(b(k,4));
    text(b(k,1)-15,b(k,2)+25,str,'FontSize',23,'color','black');
end








