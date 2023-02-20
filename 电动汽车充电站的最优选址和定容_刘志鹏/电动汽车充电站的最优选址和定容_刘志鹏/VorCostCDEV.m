function [Fcost,CCS,fcs,ucs,NchI,Ep,CVT,CVTI,CDL,CDLI]=VorCostCDEV(x,y,b,bcs,BL)
%% 说明：计算电动汽车规划目标函数，即适应值
%输入
%x，y为充电站坐标；为外界交互变换参数。x和y保证是列矩阵。
%b为充电需求点坐标和充电需求，第1,2列为横坐标和纵坐标，第3列为充电需求，为固定参数；
%bcs为配送区域的集中充电站坐标，第1,2列为横、纵坐标。
%BL为图坐标与实际坐标的比例，为固定参数。

%输出
%Fcost为年社会成本
%CCS为建设运行成本
%fcs为固定建设成本
%ucs为运行管理成本，取fcs的0.1
%NchI为各充电站规模，即充电机数量
%Ep为个充电站服务的快充车辆数，即充电需求
%CVT为用户充电途中耗时成本
%CVTI为各充电站的用户充电途中耗时成本
%CDL为电池年配送成本
%CDLI为各更换站电池配送成本


w=100*1e4;  

q=70*1e4;   

e=1.5*1e4;

beta=25;

lmd=1.2;

v=25;   

Nchmax=12;  
Nchmin=4;

ns=4;

dmax=1200;

Dmin=500;

rr=0.08;

ms=20;

Rz=(rr*(1+rr).^ms)./((1+rr).^ms-1);

mui=0.6;

dr=1;

m=1.5;


Ep=[];  


vv=VoronoiArea([x,y],3);    %V图，20次迭代，20个xy点输入,加结点10^3倍，最终由顶点生成的V图多边形结点坐标


bax=b; 
for k=1:length(x)           %x原本是个列矩阵
    Ai=find(vv(:,3)==k);
    xx=vv(Ai,1);
    yy=vv(Ai,2);
    kk=convhull(xx,yy); 
    in=inpolygon(bax(:,1),bax(:,2),xx(kk),yy(kk));
    bax(in,4)=k;
end

for i=1:length(x)
    gb=bax(bax(:,4)==i,:);
    Ep=[Ep;[sum(gb(:,3)),i]];
    
    
    dtI=sqrt(((gb(:,1))-x(i)).^2+((gb(:,2))-y(i)).^2).*BL;
    
    if all((lmd.*dtI>dmax)==0)==0
        dp=1e5;
    else
        dp=0; 
    end
   
    CVTI(i,:)=(1+dp).*365.*beta.*(sum(gb(:,3).*lmd.*dtI)./(v*1e3));
end
CVT=sum(CVTI);


NchI=round(mui.*Ep(:,1)./ns);

fcs=w+q.*NchI+e.*NchI.^2;

Np=1e3; 
NchIdex=find(NchI<Nchmin |NchI>Nchmax);

Dij=dist([x';y']).*BL;
[Da,Db]=find((lmd.*triu(Dij)>=Dmin)~=triu(ones(size(Dij)),1)); 

ULI=unique([Da;Db;NchIdex]);  

fcs(ULI)=Np.*fcs(ULI);   


ucs=0.1.*fcs;


CCS=sum(fcs.*Rz+ucs);


BcsX=bcs(:,1).*ones(length(x),1);
BcsY=bcs(:,2).*ones(length(x),1);
rij=sqrt((BcsX-x).^2+(BcsY-y).^2).*BL./1e3;
CDLI=365.*dr.*m.*Ep(:,1).*lmd.*rij;

CDL=sum(CDLI);


Fcost=CCS+CVT+CDL;