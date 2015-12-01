function dx = model_eq1_2_verA1(t,x,beta,gama,q,h,Ds,ass,k0,k1,k2,k3,k4,k5,L_mid,L_slope,z,z0)
% Implementation of the first model discussed in the paper
% `Modeling Robustness Tradeoffs in Yeast Cell Polarization Induced by
% Spatial Gradients`
% by Ching-Shan Chou, Qing Nie, Tau-Mu Yi
% please refer to equation 1.2 in the paper
global sum_a
dx = zeros(2,1);
a = x(1);
b = x(2);

sum_a = sum_a + a;
avg_a = sum_a/t;
if t < 0.001
    avg_a = 0;
end
a_kapa = avg_a-ass; 
u = L_mid+L_slope*(z-z0);

db = k4*a_kapa*b;

da_term1 = Ds*a;
da_term2 = k0/(1+(beta*u)^-q);
da_term3 = k1/(1+(gama*a)^-h);
da_term4 = k2*a;
da_term5 = k3*b*a;
da_term6 = k5*a_kapa;
da = da_term1+da_term2+da_term3-da_term4-da_term5-da_term6;

dx = [da;db];
end

