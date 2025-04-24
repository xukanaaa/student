function X=bernoulliRV(L,p)
U=rand(1,L)
X=(U<p);
end