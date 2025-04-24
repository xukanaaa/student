function [E,P]=SeqEngPower(x,ts)
E=(sum(abs(x).^2))*ts;
P=(sum(abs(x).^2))/length(x);
end