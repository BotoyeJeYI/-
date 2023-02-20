function  vv=VoronoiArea(a,en)
%a：V图顶点，第一列是x，第二列是y，行编号是区域编号；
%en: 若图边界点不在V图区域中，则伸展V图结点，增加V图区域，en表示增加结点10^en倍。
%vv：由顶点生成的V图多边形结点坐标，第1列是x，第二列y，第三列是区域编号，与a的行顺序一致。

x=a(:,1);
y=a(:,2);

[vx,vy]=voronoi(x,y);

[v,c]=voronoin([x,y]);

numInf=0; vmn1=[]; vnn=[];
for i=1:length(c)
    if any(c{i} == 1)
        numInf=numInf+1;
        iIc(numInf)=i;
        cm=c{i};
        cInf=find(cm~=1);
        cn=cm(cInf);
        vmn0=v(cn,:);
        vmn0(:,3)=i;
        vmn1=[vmn1;vmn0];
    else
        vnn0=v(c{i},:);
        vnn0(:,3)=i;
        vnn=[vnn;vnn0];
    end    
end

vxInf=vx(:,(size(vx,2)-numInf+1):size(vx,2));
vyInf=vy(:,(size(vy,2)-numInf+1):size(vy,2));

kIx=(vxInf(2,:)-vxInf(1,:))*10^en;
Ix=kIx+vxInf(1,:);
Iy=(vyInf(2,:)-vyInf(1,:))./(vxInf(2,:)-vxInf(1,:)).*(Ix-vxInf(1,:))+vyInf(1,:);

vxInf=[vxInf(1,:);Ix];
vyInf=[vyInf(1,:);Iy];


vInf=[];
for i=1:length(iIc)
    vInf0=[vxInf(2,i),vyInf(2,i),iIc(i)];
    vInf=[vInf;vInf0];
end


vre=[];
for i=1:length(iIc)
    D1=pdist([vInf(i,1:2);a(iIc(i),:)]);
    D1=round(D1*1e4)/1e4;  
    for j=i+1:length(iIc)+i-1
        if j>length(iIc)
            j=j-length(iIc);
        end
        D2=pdist([vInf(i,1:2);a(iIc(j),:)]);
        D2=round(D2*1e4)/1e4;   
        if D2==D1
            vre0=[vInf(vInf(:,3)==iIc(i),1:2),iIc(j)];
            vre=[vre;vre0];
        end
    end
end


vn1=[vmn1;vnn;vInf;vre];
vv=[];
for i=1:length(x)
    ik=find(vn1(:,3)==i);
    vv=[vv;vn1(ik,:)]; 
end