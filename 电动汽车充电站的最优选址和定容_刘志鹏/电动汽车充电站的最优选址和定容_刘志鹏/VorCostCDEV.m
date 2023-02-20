function [Fcost,CCS,fcs,ucs,NchI,Ep,CVT,CVTI,CDL,CDLI]=VorCostCDEV(x,y,b,bcs,BL)
%% ˵��������綯�����滮Ŀ�꺯��������Ӧֵ
%����
%x��yΪ���վ���ꣻΪ��罻���任������x��y��֤���о���
%bΪ������������ͳ�����󣬵�1,2��Ϊ������������꣬��3��Ϊ�������Ϊ�̶�������
%bcsΪ��������ļ��г��վ���꣬��1,2��Ϊ�ᡢ�����ꡣ
%BLΪͼ������ʵ������ı�����Ϊ�̶�������

%���
%FcostΪ�����ɱ�
%CCSΪ�������гɱ�
%fcsΪ�̶�����ɱ�
%ucsΪ���й���ɱ���ȡfcs��0.1
%NchIΪ�����վ��ģ������������
%EpΪ�����վ����Ŀ�䳵���������������
%CVTΪ�û����;�к�ʱ�ɱ�
%CVTIΪ�����վ���û����;�к�ʱ�ɱ�
%CDLΪ��������ͳɱ�
%CDLIΪ������վ������ͳɱ�


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


vv=VoronoiArea([x,y],3);    %Vͼ��20�ε�����20��xy������,�ӽ��10^3���������ɶ������ɵ�Vͼ����ν������


bax=b; 
for k=1:length(x)           %xԭ���Ǹ��о���
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