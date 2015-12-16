R = 10^-5;%sphere radius
D = 5*10^-13;%diffusion constant of the phosphorylated protein
p0 = 100;%initial concentration of protein (all unphosphorylated)
kp = 1;%rate of de-phosphorylation
ka = 1;%rate of activation (phosphorylation)

time_steps = 20;
end_time = 5;
spatial_steps = 20;