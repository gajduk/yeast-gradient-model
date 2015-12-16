defaults;
n_spatial_steps = 200;
kas = linspace(1,10000,100);

res = zeros(1,length(kas));
res_ = zeros(length(kas),n_spatial_steps);
for i=1:length(kas)
    ka = kas(i);%rate of activation- phosphorylation
    [r,t,sol,u,steady_state,half_fraction] = run_model(time_steps,end_time,n_spatial_steps,R,D,p0,ka,kp);
    res(i) = half_fraction;
    res_(i,:) = steady_state;
end

figure;
plot(r,res_(1,:));
title('ka = 1')
xlabel('Distance from nucleus [m]')
ylabel('Fraction of phophorylated protein in steady-state')
grid on

figure;
plot(r,res_(length(kas),:));
title('ka = 1000')
xlabel('Distance from nucleus [m]')
ylabel('Fraction of phophorylated protein in steady-state')
grid on

figure;
plot(kas,res/R)
xlabel('Rate of activation')
ylabel('Gradient steepness (1=very steep, 0 no gradient)')
