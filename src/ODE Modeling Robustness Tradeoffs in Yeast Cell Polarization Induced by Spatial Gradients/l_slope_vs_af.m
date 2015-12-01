%plots the concentration of a species at the shmoo pit where the ligand
%concentration is the highest for different q - ligand acivation exponents
%this is a reproduction of figure 2A in the paper
load_defaults

res = [];
qs = [1 10 100 1000];
L_slopes_exponent = -4:.1:0;
L_slopes = 10.^L_slopes_exponent;
for L_slope = L_slopes
    temp = [];
    for q = qs
        [t,y] = ode45(@(t,y) model_eq1_2_verA1(t,y,beta,gama,q,h,Ds,ass,k0,k1,k2,k3,k4,k5,L_mid,L_slope,z,z0),[0 10],y0);
        temp = [temp y(end,1)];
    end
    res = [res temp'];
end
figure;
semilogx(L_slopes,res);
xlabel('L slope')
ylabel('Species a')
legend({'q = 1','q = 10','q = 100','q = 1000'})