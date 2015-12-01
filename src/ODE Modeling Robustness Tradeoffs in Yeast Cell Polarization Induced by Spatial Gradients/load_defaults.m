global sum_a
sum_a = 0;
k0 = 1;%ligand activation of a
k1 = 0;%autocatalytic production of a
k2 = 1;%self degradation of a
k3 = 1;%integral negative feedback of b->a
k4 = 0;%autocatalytic production of b
k5 = 0;%proportional negative feeback of b->a

beta = 1;%ligand activation denominator constant
gama = 1;%autocatalytic activation denominator constant

q = 1000;%ligand activation exponent
h = 1;%autcatalytic activation exponent

ass = 1;%steady state concentration of species a
Ds = 0;%diffusion constant

R = 7;%cell radius

L_mid = 1;%ligand concentration at the middle of the cell
L_slope = 1/R;%slope of ligand concentration
z0 = 0;%coordinate of the middle of the cell

%initial conditions
a0 = 1;
b0 = 1;
y0 = [a0,b0];

z = R;%cell tip