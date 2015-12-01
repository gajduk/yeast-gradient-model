%plots the concentration of a species at different locations on the membrane
%where the ligand concentration is varies, for different beta values
%this is a reproduction of figure 2B in the paper
load_defaults
betas = [0.9,0.92,0.95,1];
res = [];
zs = -1:.05:1;
for z = zs
    temp = [];
    for beta = betas
        [t,y] = ode45(@(t,y) model_eq1_2_verA1(t,y,beta,gama,q,h,Ds,ass,k0,k1,k2,k3,k4,k5,L_mid,L_slope,z,z0),[0 10],y0);
        temp = [temp y(end,1)];
    end
    res = [res temp'];
end
figure;
plot(zs,res);
xlabel('z')
ylabel('Species a')
legend({'\beta = 0.9','\beta = 0.92','\beta = 0.95','\beta=1.0'})