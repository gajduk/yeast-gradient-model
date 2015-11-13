p = 1;%dephosphorylation rate
D = 1;%Fus3 diffusion constant
k = 1;%phosphorylation rate due to activation by ligand

init_F = 10;%initial concentration of UNphosphorylated Fus3
init_F_p = 0;%initial concentration of phosphorylated Fus3

n = 7;%num_compartments



y0 = [ones(n,1)*init_F;ones(n,1)*init_F_p];

[t,y] = ode45(@(t,y) model(t,y,p,D,k,n),[0 10],y0);

y_p = y(:,n+1:n+n)./(init_F+init_F_p);

figure;
plot(t,y_p);
xlabel('Time')
ylabel('Fraction of Fus_3  that is phosphorylated')
title(sprintf('Time profiles of %d compartment Yeast model',n))
legend_s = cell(n,1);
legend_s(1) = cellstr('Schmoo Fus3_p');
legend_s(n) = cellstr('C_L Fus3_p');
for i=2:n-1
    legend_s(i) = cellstr(sprintf('C_{%d} Fus3_p',i-1));
end
legend(legend_s);
figure;
plot(y_p(end,:))
xlabel('Distance from schmoo')
ylabel('Fraction of Fus_3 that is phosphorylated')
title(sprintf('Gradient of %d compartment Yeast model',n));
