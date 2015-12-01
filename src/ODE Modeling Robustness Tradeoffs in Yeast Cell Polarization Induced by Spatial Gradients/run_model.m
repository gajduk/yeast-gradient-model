load_defaults

[t,y] = ode45(@(t,y) model_eq1_2_verA1(t,y,beta,gama,q,h,Ds,ass,k0,k1,k2,k3,k4,k5,L_mid,L_slope,z,z0),[0 10],y0);

a = y(:,1);
b = y(:,2);

figure;
plot(t,a);
xlabel('Time')
ylabel('Species a')

figure;
plot(t,b);
xlabel('Time')
ylabel('Species b')