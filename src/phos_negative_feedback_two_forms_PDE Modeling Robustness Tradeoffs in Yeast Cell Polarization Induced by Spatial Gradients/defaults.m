R = 10^-5;%sphere radius
D = 10*10^-13;%diffusion constant of the phosphorylated protein
p0 = 300;%[nM]initial concentration of protein (all unphosphorylated)
kp = 1;%s^-1 rate of de-phosphorylation
ka = 1;%s^-1 rate of activation (phosphorylation)

%phosphatase governed by Michaelis Menten kinetics
k_cat = .05;%s^-1 catalytic activity of the phosphatases
E0 = 75;%[nM] is the initial concentration of the phosphotases
K_M = 100;%[nM] is the  Michaelis constant of the phosphatases
A = 10;%strength of the feedback (a free parameter).

time_steps = 50;
end_time = 30;
spatial_steps = 40;