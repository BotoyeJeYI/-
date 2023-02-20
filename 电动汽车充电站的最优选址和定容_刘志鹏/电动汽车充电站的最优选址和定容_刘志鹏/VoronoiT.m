function [vxT,vyT]=VoronoiT(x,y,en)

mx=x(1)+(x(2)-x(1))./2; 
my=y(1)+(y(2)-y(1))./2; 
mk=-(x(2)-x(1))./(y(2)-y(1));
vxT=[x(1)-(x(2)-x(1)).*en
    x(2)+(x(2)-x(1)).*en];
vyT=mk.*(vxT-mx)+my;