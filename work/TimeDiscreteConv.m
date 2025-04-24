function [x,n]=TimeDiscreteConv(x1,n1,x2,n2)
x=conv(x1,x2);
k0=n1(1)+n2(1);
L=length(x1)+length(x2)-1;
n=k0:k0+L-1;
end