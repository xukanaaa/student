function [t,y]=sampling(f,C,Ts)
t=0:Ts:C-Ts;
y=eval(subs(f,'t',t));
end