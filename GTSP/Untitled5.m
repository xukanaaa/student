phySkew=80e-6;%物理时钟频偏最大值
alpha1=1-phySkew+(2*phySkew)*rand;
alpha2=1-phySkew+(2*phySkew)*rand;
arr=zeros(10000,1);
for i=1:10000
    arr(i)=(alpha1+2e-6*randn+exprnd(4e-6)+2e-6*randn+exprnd(4e-6))/(alpha2+2e-6*randn+exprnd(4e-6)+2e-6*randn+exprnd(4e-6));
end
alpha=alpha1/alpha2;
w=0;
for i=10:10000
    Alpha=w*(arr(i-9)+arr(i-8)+arr(i-7)+arr(i-6)+arr(i-5)+arr(i-4)+arr(i-3)+arr(i-2)+arr(i-1)+arr(i));
    w=w+0.01*(alpha-Alpha)*arr(i);
end
disp(w);
