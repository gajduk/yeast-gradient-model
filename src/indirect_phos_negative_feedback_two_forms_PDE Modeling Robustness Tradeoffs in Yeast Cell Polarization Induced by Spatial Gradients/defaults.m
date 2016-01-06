R = 10^-5;%sphere radius
D = 5*10^-14;%diffusion constant of the phosphorylated protein
p0 = 300;%[nM]initial concentration of protein (all unphosphorylated)
phos0 = 50;%[nM] initial concentration of phospatases
phos_max = phos0*4;%[nM] maximal concentration of phosphatases
phos_regulation_rate = .002;% - the speed of phosphatase dynamics

kp = .003;%s^-1 rate of de-phosphorylation
ka = .008;%s^-1 rate of activation (phosphorylation)

A = 1;% strength of the feedback (a free parameter).


time_steps = 3000;
end_time = 3000;
spatial_steps = 40;