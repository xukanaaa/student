alpha_i=1.00005;
alpha_j=0.99994;
beta_i=0.78;
beta_j=-0.5;
round=100;
T_send=zeros(round,1);
T_receive=zeros(round,1);
alpha_ij=zeros(round,1);
alpha_ij(1,1)=alpha_i/alpha_j;
delta_alpha=zeros(round,1);

for i=1:round
    T_send(i,1)=alpha_i*i+beta_i+5e-9*randn;
    T_receive(i,1)=alpha_j*i+beta_j+300/3e8+5e-9*randn;
    if i>=2
    alpha_ij(i,1)=(T_send(i,1)-T_send(1,1))/(T_receive(i,1)-T_receive(1,1));
    delta_alpha(i,1)=(alpha_ij(i,1)-alpha_ij(1,1))*(i-1);
    end
end

format long;
disp(delta_alpha);


