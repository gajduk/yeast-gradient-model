function [r,t,sol,u,steady_state,half_fraction] = run_model(time_steps,end_time,spatial_steps,R,D,p0,ka,kp,k_cat,E0,K_M,A)
    %for a description of the model see readme
    %
    %parameters
    %time_steps,spatial_steps - number of discretization steps in time,space respectively (higher values better precision and
    %smoother solutions but longer running time
    %end_time - the system solutions at times 0 [s] - end_time [s] are
    %calculated
    %R - cell radius
    %D - diffusion constant of the phosphorylated protein
    %p0 - initial concentration of protein (all unphosphorylated)
    %kp - rate of de-phosphorylation
    %ka - rate of activation (phosphorylation)
    %
    %phosphatases Michaelis Menten kinetics
    %k_cat  catalytic activity of the phosphatases
    %E_0  is the initial concentration of the phosphotases
    %K_M  is the  Michaelis constant of the phosphatases
    %A is the strength of the feedback (a free parameter).
    
    %returns
    %r - the spatial grid (one dimensional array with `spatial_steps` elements - elements E [0,R]
    %t - time grid (one dimensional array with `time_steps` elements- elements E [0,end_time]
    %sol - the solution of the system of PDEs for each point - (r,t) (three
    %dimensional matrix - where axis1 = space, axis2 = time and axis3 =
    %variable (unphosphorylated or phosphorylated protein)
    %u - fraction of phosphorulated protein at each point (r,t) (two
    %dimensional matrix with values E [0,1]
    %steady state - the fraction of phosphorylated proteins at t=end_time
    %(one dimensional array - each value corresponds to a point in space)
    %half_fraction - the value of `r`/R where the fraction of phosphorylated
    %protein is half of the fraction at the cell membrane e [0,1] - values
    %of 1 would mean an extreme gradient where al the phosphorylated
    %protein is at the membrane, and value of 0 would mean that
    %phosphorylated protein is uniformly distributed throught the cell
        
    xmesh = linspace(0,R,spatial_steps);
    tspan = linspace(0,end_time,time_steps);
    
    m = 2;
    pdefun = @(x,t,u,DuDx) pdex1pde(x,t,u,DuDx,D,kp,ka,R,k_cat,E0,K_M,A);%actual PDE
    icfun = @(x) pdex1ic(x,p0);%initial conditions
    bcfun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t);%boundary conditions

    sol = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan);

    u = sol(:,:,2)./(sol(:,:,2)+sol(:,:,1));


    steady_state = u(end,:);
    half = 1;
    for k=1:length(xmesh)
        if steady_state(k) < steady_state(end)/2
            half = k;
        end
    end
    half_fraction = xmesh(half)/R;
    r = xmesh;
    t = tspan;
end