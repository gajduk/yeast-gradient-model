
R = 10^-5;%sphere radius
D = 1*10^-13;%diffusion constant of the phosphorylated protein
pR = 1;%initial condition of the phosphorylated protein at the membrane p(x=R,t=0)
kp = 0;%rate of de-phosphorylation

xmesh = linspace(0,R,100);
tspan = linspace(0,100,10);

m = 2;

pdefun = @(x,t,u,DuDx) pdex1pde(x,t,u,DuDx,D,kp);%actual PDE
icfun = @(x) pdex1ic(x);%initial conditions
bcfun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,D,pR);%boundary conditions

sol = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan);

u = sol(:,:,1);
    
surf(xmesh,tspan,u)    
title('Numerical solution')
xlabel('Distance x [m]')
ylabel('Time t [s]')
zlabel('Concentration of phosphorylated protein')