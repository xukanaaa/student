x=zeros(1,50);
for i=1:50
    x(1,i)=(5*randn+20)*10^(-6)*((-1)^randi([0,1]));
end
disp(x);