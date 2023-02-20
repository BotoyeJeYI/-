clear ; 
clc; 
close all
%% ��������
%    1    2     3         4        5
%b=[�������꣬���󸺺ɣ��ɱ�ϵ��������]
%������������
b(:,1) =1.0e+003.*([0.1092, 0.1197,0.2578,0.4259, 0.1257,0.2803,0.4439,0.5505,0.5610,0.5700,0.1332,0.3013,0.4559,0.5850,0.1452,0.3163,0.4739,0.5880,0.3193,0.4784,0.6015,0.6736,0.8657,1.0308,0.6781,0.8327,1.0188,0.6811,0.8192,1.0083,0.6691,0.6976,0.8191,1.0098]');
b (:,2)=1.0e+003.*[0.1348;0.2399;0.1724;0.1739;0.3420;0.3375;0.3360;0.1108;0.2024;0.3075;0.4591;0.4455;0.4380;0.4260;0.6092;0.5341;0.5341;0.5341;0.6407;0.6452; 0.6467;0.1574 ;0.1649;0.1634;0.2924;0.2909;0.2939;0.4425; 0.4576;0.4546;0.5431;0.6392;0.6377;0.6347];
%�������㳣��������ɵ㸺��b(:,3)��kW��
b(:,3)=[2480;2480;8680;11400;890;2340;4160;560;1670;5010;2670;8280;7400;1430;7500;4840;3400;4290;3840;3680;2560;7000;14800;8960;3160;7000;5000;2280;10360;10000;760;6000;7040;5600];

%���г��վ����,������������һ��
bcs=[ 937.7296  379.5010;
  310.3141  238.4076];
%ѡַ����
Tn=7;   

%% �ɱ�����
na=4500;    %�綯��������
alp=0.1;  
b(:,4)=round(alp.*b(:,3)./sum(b(:,3)).*na); %ÿ�������ƽ������
b(23,4)=37;  
ns=4;
mui=0.6;   
Nchz=round(mui.*sum(b(:,4))./ns);
bm=1.0e+003*[0.0086,0.0088;1.1734,0.0088;1.1734,0.7412;0.0086,0.7412;0.0086,0.0088];
BL=sqrt(8.2*1.0e6./((max(bm(:,1))-min(bm(:,1)))*(max(bm(:,2))-min(bm(:,2)))));%BLΪͼ������ʵ������ı�����Ϊ�̶�����

%% ���򻮷�
Area1=1.0e+003 *[0.9377   -1.0860;
1.1734    0.0088;
1.1734    0.7412;
0.3103    1.7040;
0.9377   -1.0860];
Area1=[Area1,1.*ones(size(Area1,1),1)]; %����size(Area1,1)�������Area1������
Area2=1.0e+003 *[0.0086    0.0088;
0.9377   -1.0860;
0.3103    1.7040;
0.0086    0.7412;
0.0086    0.0088];
Area2=[Area2,2.*ones(size(Area2,1),1)];
%����ֽ�
vv=[Area1;Area2];   %10*3����
for k=1:size(bcs,1) %size(bcs,1)���bcs��������2,��k=1:2
    Ai=find(vv(:,3)==k);    %��vv�ĵ����в��ҵ���k��Ԫ�ط�������ֵ,���б�vvԪ���������򼸣�����Ϊ������ʱ����Ai��
    xx=vv(Ai,1);            %�����꣬�������㸺�ɣ�����k�ĵ� ������
    yy=vv(Ai,2);            %������
    kk=convhull(xx,yy);     %����͹����kk��һ��������
    %in = inpolygon(x,y,xv,yv)%ע��xv,yv�����˶���α߽硣x,y��Ӧ���ǵ������꣬�ж��Ƿ��ڶ�����ڡ�
    %���ؽ��Ϊ�߼�logical���ͣ�������������Ŷ��������ڶ�Ӧ�ľͷ���1������Ϊ0��
    in=inpolygon(b(:,1),b(:,2),xx(kk),yy(kk));
    b(in,5)=k; 
end
%���������1���ɱ�ϵ��֮�ͣ�2��mui*�ɱ�ϵ��֮��/ns��ֵȡ����3������
Ep=[];
for i=1:size(bcs,1)     %2�����򣬼�i=1��2
    gb=b(b(:,5)==i,:);  %�������ڵ������ݴ浽gb
    Ep=[Ep;[sum(gb(:,4)),round(mui.*sum(gb(:,4))./ns),i]];  %������������ݼ��㼰����
end

%% ����Ⱥ�㷨����

PopSize=20; %��Ⱥ����
MaxIter=300;  %����������
%ѧϰ����
c1s=2.5;    %����ѧϰ�������ֵ
c2s=0.5;    %���ѧϰ�������ֵ
c1e=0.5;    %����ѧϰ������Сֵ
c2e=2.5;    %���ѧϰ������Сֵ
w_start=0.9;   %Ȩֵ���ֵ
w_end=0.4;     %Ȩֵ��Сֵ
Iter=1;        

%������Լ��������������
xmax=max(Area1(:,1));   %����1���
xmin=min(Area1(:,1));   %����1��С��
ymax=max(Area1(:,2));   %����2����
ymin=min(Area1(:,2));   %����2��С��
x = xmin + (xmax-xmin).*rand(Tn,PopSize);   %����x��������ӣ�Tn��ѡַ+PopSize����Ⱥ
y = ymin + (ymax-ymin).*rand(Tn,PopSize);   %����y���������
%% ��ʼ��
X=[x;y];    %�ϲ�xy���ӣ��γ�һ�������ĸ���,ǰ7�б�ʾX���꣬��7�б�ʾY���꣬��������������ʾ��������ʹ��
V=rand(Tn*2,PopSize);   %������������ٶ�
Vmax=0.1*max((xmax-xmin),(ymax-ymin));      %��������ٶ�
inAr1=find(b(:,5)==1);                      %�ҵ�����1�ı�ţ����������洢��inAr1��
bb=[b(inAr1,1:2),b(inAr1,4)];               %bb�ݴ���������ͳɱ�ϵ��
%��ʼѰ��
for pk=1:1:PopSize
    [FX(pk),~,~,~,~,~,~,~,~,~]=VorCostCDEV(X(1:Tn,pk),X(Tn+1:end,pk),bb,bcs(1,:),BL);   %����Ӧ��ֵ,~�������ĳ������������ֻ�����һ��������X(1:Tn,pk)��ʾpk����Ⱥ��Tn������ѡַ�ĺ�����,X(Tn+1:end,pk)��ʾpk����Ⱥ��Tn������ѡַ��������
end

%�����ʼ����
PBest=X;            %PBest��ʼ��������
FPBest=FX;          %FPBest��ʼ����Ŀ�꺯��ֵ


[Fgbest,r]=min(FX); %����Ŀ�꺯��ֵ�ĳ�ʼ��Сֵ����Сֵ���ڼ��У����ڼ�����Ⱥ����
Fm(Iter)=Fgbest;    %������Ŀ�꺯��ֵ�ĳ�ʼ��Сֵ����Fm��
CF=Fgbest;          %������Ŀ�꺯��ֵ�ĳ�ʼ��Сֵ����CF��
Best=X(:,r);        %��ʼ������Ⱥ���������Best����������ǰ7�����꣬��7�����꣩
FBest(Iter)=Fgbest; %������Ŀ�꺯��ֵ�ĳ�ʼ��Сֵ����FBest��1����

FgNum=0;
%% ����ȺѰ��
while (Iter<=MaxIter)   %����Ⱥ�㷨
    Iter=Iter+1;        %��������
    w_now=((w_start-w_end)*(MaxIter-Iter)/MaxIter)+w_end;   %����Ȩ�� 
    A=repmat(X(:,r),1,PopSize);     
    R1=rand(Tn*2,PopSize);
    R2=rand(Tn*2,PopSize);   
    c1=c1e+(c1s-c1e)*(1-acos(-2*Iter/(MaxIter+1)+1)/pi);
    c2=c2e+(c2s-c2e)*(1-acos(-2*Iter/(MaxIter+1)+1)/pi);    
    V=w_now*V+c1*R1.*(PBest-X)+c2*R2.*(A-X);                %�����ٶȸ��¹�ʽ
    changeRows=V>Vmax;
    V(changeRows)=Vmax;   
    changeRows=V<-Vmax;
    V(changeRows)=-Vmax;    
    X=X+1.0*V;  
   for pk=1:1:PopSize
       [FX(pk),~,~,~,~,~,~,~,~,~]=VorCostCDEV(X(1:Tn,pk),X(Tn+1:end,pk),bb,bcs(1,:),BL);    %����Ӧ��
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

%% �������ֵ
 %����
x =Best(1:Tn);
y =Best(Tn+1:end);
[f1,f2,f3,f4,Num_chI,num_Epc,f5,f6,f7,f8]=VorCostCDEV(x,y,bb,bcs(1,:),BL)
FBest  
Best 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��ͼ
%% ����ͼ
figure(1)
Iter_t=1:1:MaxIter+1;
plot(Iter_t,Fm,'.-')
legend('ÿ�ε���ֵ')
title('����ͼ')

%% ѡַͼ
figure(2)
a=imread('T3.Png');
imshow(a);   
hold on
[vxT,vyT]=VoronoiT(bcs(:,1),bcs(:,2),0);  
%����ʽ���վ
plot(bcs(:,1),bcs(:,2),'gs','linewidth',12);
%����ֽ�
plot(vxT,vyT,'r -','linewidth',3);   
%�����
plot(b(:,1),b(:,2),'k*','linewidth',5)  ;
%�ܱ�
plot(bm(:,1),bm(:,2),'k-','linewidth',3) 
%���������
for k=1:length(b(:,4))
    str=num2str(b(k,4));
    text(b(k,1)-15,b(k,2)+25,str,'FontSize',23,'color','black');
end
axis equal
 [vx,vy]=voronoi(x,y);
 plot(x,y,'b^','linewidth',3); 
 %���ų��վ����
for k=1:length(x)
    str=num2str(k);
    text(x(k),y(k),str,'FontSize',20,'color','red');
 end
axis([min(bm(:,1))-3 max(bm(:,1))+3 min(bm(:,2))-3 max(bm(:,2))+3])
legend('���г��վ','����ֽ�','��������','�ܱ߽���','���վ����ѡַ��')
title('�滮ѡַͼ')

%% ѡַͼ
figure(3)
a=imread('T3.bmp');
imshow(a);   
hold on
[vxT,vyT]=VoronoiT(bcs(:,1),bcs(:,2),0);  
%����ʽ���վ
plot(bcs(:,1),bcs(:,2),'gs','linewidth',12);
%����ֽ�
plot(vxT,vyT,'r -','linewidth',3);   
%�����
plot(b(:,1),b(:,2),'k*','linewidth',5)  ;
%�ܱ�
plot(bm(:,1),bm(:,2),'k-','linewidth',3) 
%���������
for k=1:length(b(:,4))
    str=num2str(b(k,4));
    text(b(k,1)-15,b(k,2)+25,str,'FontSize',23,'color','black');
end








